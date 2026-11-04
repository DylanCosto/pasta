//! The normalization of parsed SAM records that HTSlib's `sam_parse1` does.
//!
//! HTSlib makes the placement fields of a SAM record consistent as it parses it, and so does
//! pasta:
//!
//! - An RNAME not in the header (by `SN` or an `AN` alternative name) is treated as `*`, with a
//!   warning. An RNAME other than `*` is an error if the header has no `@SQ` lines.
//! - An RNAME with POS 0 is treated as `*`, with a warning.
//! - A record with no reference sequence gets the unmapped flag (0x4).
//! - A record with CIGAR `*` gets the unmapped flag, with a warning if it lacked it.
//! - An RNEXT not in the header is treated as `*`, with a warning.
//! - An RNEXT with PNEXT 0 is treated as `*`, with a warning.
//!
//! Where HTSlib's result is wrong, pasta differs:
//!
//! - A position whose reference sequence is dropped because its name is not in the header is
//!   dropped too: POS (for RNAME) and PNEXT (for RNEXT, or RNEXT `=` of such an RNAME) become 0.
//!   HTSlib keeps them, giving records like `* 10` that place a read on no reference sequence
//!   (and that htsjdk rejects: "Alignment start should be 0 because reference name = *").
//! - RNEXT `=` means the RNAME as written. HTSlib resolves it after dropping an RNAME with POS 0,
//!   and so also drops the mate's reference sequence (leaving `* PNEXT`); pasta keeps the mate
//!   on that reference sequence.
//! - HTSlib warns for every record; each distinct warning is printed once per reader.

use std::{collections::HashMap, io};

use bstr::{BStr, BString, ByteSlice};

use crate::sam::{
    Header,
    alignment::{RecordBuf, record::Flags},
    header::record::value::map::reference_sequence::tag::ALTERNATIVE_NAMES,
    record::Fields,
};
use crate::variant::log;

const MISSING: &[u8] = b"*";
const EQ: &[u8] = b"=";

// Mandatory field indices.
const FLAG: usize = 1;
const RNAME: usize = 2;
const POS: usize = 3;
const RNEXT: usize = 6;
const PNEXT: usize = 7;

/// A zero (unavailable) POS or PNEXT.
const ZERO: &[u8] = b"0";

/// The function HTSlib attributes these messages to.
const FUNC: &str = "sam_parse1";

/// The `hts_strprint` buffer `sam_parse1` quotes reference sequence names into.
const NAME_BUF_LEN: usize = 40;

/// The state HTSlib's record normalization needs across records: the header's alternative
/// reference sequence names and the warnings already printed.
#[derive(Debug, Default)]
pub(crate) struct Normalizer {
    /// Alternative name (`@SQ AN`) to reference sequence name, built from the header on the first
    /// name that is not a reference sequence name.
    alternative_names: Option<HashMap<BString, BString>>,
    warned_names: HashMap<BString, Warned>,
    warned: Warned,
}

#[derive(Clone, Copy, Debug, Default)]
struct Warned {
    reference: bool,
    mate: bool,
    position: bool,
    cigar: bool,
    mate_position: bool,
}

/// A reference sequence field as resolved by `bam_name2id`.
enum Resolved<'a> {
    Missing,
    /// The ID and the name to keep: `None` if the field is the reference sequence name itself,
    /// otherwise the reference sequence name of an alternative name.
    Found(usize, Option<&'a BStr>),
    Unknown,
}

impl Normalizer {
    /// Forgets the header, e.g., after another header is read.
    pub(crate) fn reset(&mut self) {
        *self = Self::default();
    }

    fn resolve<'h>(&mut self, header: &'h Header, name: &[u8]) -> Resolved<'h> {
        let reference_sequences = header.reference_sequences();

        if let Some(id) = reference_sequences.get_index_of(name) {
            return Resolved::Found(id, None);
        }

        let alternative_names = self
            .alternative_names
            .get_or_insert_with(|| build_alternative_names(header));

        alternative_names
            .get(name.as_bstr())
            .and_then(|canonical| reference_sequences.get_full(canonical.as_slice()))
            .map_or(Resolved::Unknown, |(id, canonical, _)| {
                Resolved::Found(id, Some(canonical.as_bstr()))
            })
    }

    fn resolve_field<'h>(
        &mut self,
        header: &'h Header,
        src: &[u8],
        is_mate: bool,
    ) -> io::Result<Resolved<'h>> {
        if src == MISSING {
            return Ok(Resolved::Missing);
        }

        if !is_mate && header.reference_sequences().is_empty() {
            log::error(FUNC, format_args!("no SQ lines present in the header"));

            return Err(io::Error::new(
                io::ErrorKind::InvalidData,
                "no SQ lines present in the header",
            ));
        }

        let resolved = self.resolve(header, src);

        if let Resolved::Unknown = resolved {
            self.warn_unknown_name(src, is_mate);
        }

        Ok(resolved)
    }

    fn warn_unknown_name(&mut self, name: &[u8], is_mate: bool) {
        let warned = self.warned_names.entry(name.into()).or_default();
        let flag = if is_mate {
            &mut warned.mate
        } else {
            &mut warned.reference
        };

        if !*flag {
            *flag = true;

            let kind = if is_mate {
                "mate reference"
            } else {
                "reference"
            };

            log::warn(
                FUNC,
                format_args!(
                    "unrecognized {kind} name {}; treated as unmapped",
                    log::strprint(name, b'"', NAME_BUF_LEN)
                ),
            );
        }
    }

    fn warn(&mut self, select: fn(&mut Warned) -> &mut bool, message: &str) {
        let flag = select(&mut self.warned);

        if !*flag {
            *flag = true;
            log::warn(FUNC, format_args!("{message}"));
        }
    }

    fn warn_position(&mut self) {
        self.warn(
            |w| &mut w.position,
            "mapped query cannot have zero coordinate; treated as unmapped",
        );
    }

    fn warn_cigar(&mut self) {
        self.warn(
            |w| &mut w.cigar,
            "mapped query must have a CIGAR; treated as unmapped",
        );
    }

    fn warn_mate_position(&mut self) {
        self.warn(
            |w| &mut w.mate_position,
            "mapped mate cannot have zero coordinate; treated as unmapped",
        );
    }

    /// Normalizes the placement fields of a lazily parsed SAM record, rewriting its text where
    /// they change.
    pub(crate) fn normalize_fields(
        &mut self,
        header: &Header,
        fields: &mut Fields,
    ) -> io::Result<()> {
        let flags = Flags::from(fields.flags()?);
        let bounds = &fields.bounds;
        let src = &fields.buf[bounds.reference_sequence_name_range()];

        // RNAME and POS. `name` is what RNEXT `=` refers to: the reference sequence name (an
        // alternative name's) if the RNAME is in the header.
        let resolved = self.resolve_field(header, src, false)?;

        let (name, mut reference_sequence_name, is_unknown) = match resolved {
            Resolved::Missing => (None, None, false),
            Resolved::Found(_, None) => (Some(src.to_vec()), None, false),
            Resolved::Found(_, Some(canonical)) => {
                (Some(canonical.to_vec()), Some(canonical.to_vec()), false)
            }
            Resolved::Unknown => (None, Some(MISSING.to_vec()), true),
        };

        let mut is_placed = name.is_some();

        if is_placed && is_zero(&fields.buf[bounds.alignment_start_range()]) {
            self.warn_position();
            is_placed = false;
            reference_sequence_name = Some(MISSING.to_vec());
        }

        let mut new_flags = flags;

        if !is_placed {
            new_flags.insert(Flags::UNMAPPED);
        }

        // CIGAR
        if fields.buf[bounds.cigar_range()].starts_with(MISSING) {
            if !new_flags.is_unmapped() {
                self.warn_cigar();
            }

            new_flags.insert(Flags::UNMAPPED);
        }

        // RNEXT and PNEXT
        let src = &fields.buf[bounds.mate_reference_sequence_name_range()];

        let (is_mate_placed, mut mate_reference_sequence_name, is_mate_unknown) = if src == EQ {
            match name {
                // `=` stays `=` while RNAME stays as written.
                Some(_) if is_placed => (true, None, false),
                Some(name) => (true, Some(name), false),
                None if is_unknown => (false, Some(MISSING.to_vec()), true),
                None => (false, Some(MISSING.to_vec()), false),
            }
        } else {
            match self.resolve_field(header, src, true)? {
                Resolved::Missing => (false, None, false),
                Resolved::Found(_, name) => (true, name.map(|name| name.to_vec()), false),
                Resolved::Unknown => (false, Some(MISSING.to_vec()), true),
            }
        };

        if is_mate_placed && is_zero(&fields.buf[bounds.mate_alignment_start_range()]) {
            self.warn_mate_position();
            mate_reference_sequence_name = Some(MISSING.to_vec());
        }

        // Rewrite from the last field changed to the first, so that earlier ranges stay valid.
        if is_mate_unknown {
            fields.replace(PNEXT, ZERO);
        }

        if let Some(name) = mate_reference_sequence_name {
            fields.replace(RNEXT, &name);
        }

        if is_unknown {
            fields.replace(POS, ZERO);
        }

        if let Some(name) = reference_sequence_name {
            fields.replace(RNAME, &name);
        }

        if new_flags != flags {
            let buf = u16::from(new_flags).to_string();
            fields.replace(FLAG, buf.as_bytes());
        }

        Ok(())
    }

    /// Resolves a reference sequence name field of a record parsed into a buffer, returning the
    /// ID and whether the name was dropped because it is not in the header.
    pub(crate) fn reference_sequence_id(
        &mut self,
        header: &Header,
        src: &[u8],
        is_mate: bool,
    ) -> io::Result<(Option<usize>, bool)> {
        Ok(match self.resolve_field(header, src, is_mate)? {
            Resolved::Found(id, _) => (Some(id), false),
            Resolved::Missing => (None, false),
            Resolved::Unknown => (None, true),
        })
    }

    /// Normalizes the placement fields of a record parsed into a buffer.
    pub(crate) fn normalize_record_buf(&mut self, record: &mut RecordBuf, is_cigar_missing: bool) {
        if record.reference_sequence_id().is_some() && record.alignment_start().is_none() {
            self.warn_position();
            *record.reference_sequence_id_mut() = None;
        }

        if record.reference_sequence_id().is_none() {
            record.flags_mut().insert(Flags::UNMAPPED);
        }

        if is_cigar_missing {
            if !record.flags().is_unmapped() {
                self.warn_cigar();
            }

            record.flags_mut().insert(Flags::UNMAPPED);
        }
    }

    /// Normalizes the mate fields of a record parsed into a buffer.
    pub(crate) fn normalize_record_buf_mate(&mut self, record: &mut RecordBuf) {
        if record.mate_reference_sequence_id().is_some() && record.mate_alignment_start().is_none()
        {
            self.warn_mate_position();
            *record.mate_reference_sequence_id_mut() = None;
        }
    }
}

/// Returns whether a POS or PNEXT field is 0, which HTSlib reads as unplaced. Only `0` itself
/// counts, as pasta parses positions (`00` is an invalid position, not 0).
fn is_zero(src: &[u8]) -> bool {
    src == ZERO
}

/// Maps each `@SQ AN` alternative name to its reference sequence name, as HTSlib's
/// `sam_hrecs_add_ref_altnames` adds them: comma-separated, empty names skipped, the first
/// reference sequence to claim a name keeping it, and reference sequence names taking precedence.
fn build_alternative_names(header: &Header) -> HashMap<BString, BString> {
    let mut alternative_names = HashMap::new();

    for (name, reference_sequence) in header.reference_sequences() {
        let Some(list) = reference_sequence.other_fields().get(&ALTERNATIVE_NAMES) else {
            continue;
        };

        for alternative_name in list.split(|&b| b == b',').filter(|s| !s.is_empty()) {
            alternative_names
                .entry(alternative_name.into())
                .or_insert_with(|| name.clone());
        }
    }

    alternative_names
}
