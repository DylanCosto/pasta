//! SAM header.
//!
//! A SAM header is a list of header records. There are four map-like record types ([header]
//! (`HD`), [reference sequence] (`SQ`), [read group] (`RG`), and [program] (`PG`)) and one string
//! type (comment (`CO`)).
//!
//! Records are optional, which means an empty header is considered a valid SAM header.
//!
//! If there is a header (`HD`) record, it must appear as the first record.
//!
//! Reference sequence, read group, and program records are maps of records with unique
//! identifiers. Comments are a list of strings. Reference sequences must be ordered; whereas read
//! groups, programs, and comments can be unordered. (`sam::Header` defines them to be ordered.)
//!
//! [header]: `record::value::map::Header`
//! [reference sequence]: `ReferenceSequence`
//! [read group]: `ReadGroup`
//! [program]: `record::value::map::Program`
//!
//! # Examples
//!
//! ## Parse a SAM header
//!
//! ```
//! use pasta::sam;
//!
//! let s = "\
//! @HD\tVN:1.6\tSO:coordinate
//! @SQ\tSN:sq0\tLN:8
//! @SQ\tSN:sq1\tLN:13
//! ";
//!
//! let header: sam::Header = s.parse()?;
//!
//! assert!(header.header().is_some());
//! assert_eq!(header.reference_sequences().len(), 2);
//! assert!(header.read_groups().is_empty());
//! assert!(header.programs().as_ref().is_empty());
//! assert!(header.comments().is_empty());
//! # Ok::<(), sam::header::ParseError>(())
//! ```
//!
//! ## Create a SAM header
//!
//! ```
//! use std::num::NonZero;
//!
//! use pasta::sam::{
//!     self as sam,
//!     header::record::value::{map::ReferenceSequence, Map},
//! };
//!
//! let header = sam::Header::builder()
//!     .set_header(Default::default())
//!     .add_reference_sequence(
//!         "sq0",
//!         Map::<ReferenceSequence>::new(NonZero::try_from(8)?),
//!     )
//!     .add_reference_sequence(
//!         "sq1",
//!         Map::<ReferenceSequence>::new(NonZero::try_from(13)?),
//!     )
//!     .build();
//!
//! assert!(header.header().is_some());
//! assert_eq!(header.reference_sequences().len(), 2);
//! assert!(header.read_groups().is_empty());
//! assert!(header.programs().as_ref().is_empty());
//! assert!(header.comments().is_empty());
//! # Ok::<(), Box<dyn std::error::Error>>(())
//! ```

mod builder;
pub(crate) mod layout;
mod parser;
mod programs;
pub mod record;

pub use self::{
    builder::Builder,
    parser::{ParseError, Parser},
    record::Record,
};

use std::{
    io,
    str::{self, FromStr},
    sync::Arc,
};

use bstr::BString;
use indexmap::IndexMap;

use self::layout::Layout;
pub use self::programs::Programs;
use self::record::value::{
    Map,
    map::{self, Program, ReadGroup, ReferenceSequence},
};

/// A reference sequence dictionary.
pub type ReferenceSequences = IndexMap<BString, Map<ReferenceSequence>>;

/// An ordered map of read groups.
pub type ReadGroups = IndexMap<BString, Map<ReadGroup>>;

/// A SAM header.
///
/// Records are grouped by their types: header, reference sequence, read group, program, and
/// comment.
///
/// A header that was read also keeps its lines as read, which writing it reproduces, as htslib
/// does: lines in their order, tags as written, and the lines that the records do not hold (a
/// later `@HD`, a later `@RG` or an earlier `@PG` with the same ID as another, and lines of
/// other types). Records that change are rewritten in place, and new records follow the last
/// line of their type. These lines are not part of the header's value: two headers with the
/// same records are equal.
#[derive(Clone, Debug, Default)]
pub struct Header {
    header: Option<Map<map::Header>>,
    reference_sequences: ReferenceSequences,
    read_groups: ReadGroups,
    programs: Programs,
    comments: Vec<BString>,
    layout: Option<Arc<Layout>>,
}

impl PartialEq for Header {
    fn eq(&self, other: &Self) -> bool {
        self.header == other.header
            && self.reference_sequences == other.reference_sequences
            && self.read_groups == other.read_groups
            && self.programs == other.programs
            && self.comments == other.comments
    }
}

impl Eq for Header {}

impl Header {
    /// Returns a builder to create a SAM header.
    ///
    /// # Examples
    ///
    /// ```
    /// use pasta::sam;
    /// let builder = sam::Header::builder();
    /// ```
    pub fn builder() -> Builder {
        Builder::default()
    }

    /// Returns the SAM header header if it is set.
    ///
    /// # Examples
    ///
    /// ```
    /// use pasta::sam::{
    ///     self as sam,
    ///     header::record::value::{map::Header, Map},
    /// };
    ///
    /// let header = sam::Header::default();
    /// assert!(header.header().is_none());
    ///
    /// let header = sam::Header::builder()
    ///     .set_header(Map::<Header>::default())
    ///     .build();
    ///
    /// assert!(header.header().is_some());
    /// ```
    pub fn header(&self) -> Option<&Map<map::Header>> {
        self.header.as_ref()
    }

    /// Returns a mutable reference to the SAM header header if it is set.
    ///
    /// # Examples
    ///
    /// ```
    /// use pasta::sam::{self as sam, header::record::value::{map, Map}};
    ///
    /// let mut header = sam::Header::default();
    /// assert!(header.header().is_none());
    ///
    /// let hd = Map::<map::Header>::default();
    /// *header.header_mut() = Some(hd.clone());
    /// assert_eq!(header.header(), Some(&hd));
    /// ```
    pub fn header_mut(&mut self) -> &mut Option<Map<map::Header>> {
        &mut self.header
    }

    /// Returns the SAM header reference sequences.
    ///
    /// This is also called the reference sequence dictionary.
    ///
    /// # Examples
    ///
    /// ```
    /// use std::num::NonZero;
    ///
    /// use pasta::sam::{
    ///     self as sam,
    ///     header::record::value::{map::ReferenceSequence, Map},
    /// };
    ///
    /// let header = sam::Header::builder()
    ///     .add_reference_sequence(
    ///         "sq0",
    ///         Map::<ReferenceSequence>::new(NonZero::try_from(13)?)
    ///     )
    ///     .build();
    ///
    /// let reference_sequences = header.reference_sequences();
    /// assert_eq!(reference_sequences.len(), 1);
    /// assert!(reference_sequences.contains_key(&b"sq0"[..]));
    /// # Ok::<(), Box<dyn std::error::Error>>(())
    /// ```
    pub fn reference_sequences(&self) -> &ReferenceSequences {
        &self.reference_sequences
    }

    /// Returns a mutable reference to the SAM header reference sequences.
    ///
    /// This is also called the reference sequence dictionary.
    ///
    /// # Examples
    ///
    /// ```
    /// use std::num::NonZero;
    ///
    /// use pasta::sam::{
    ///     self as sam,
    ///     header::record::value::{map::ReferenceSequence, Map},
    /// };
    ///
    /// let mut header = sam::Header::default();
    ///
    /// header.reference_sequences_mut().insert(
    ///     String::from("sq0").into(),
    ///     Map::<ReferenceSequence>::new(NonZero::try_from(13)?),
    /// );
    ///
    /// let reference_sequences = header.reference_sequences();
    /// assert_eq!(reference_sequences.len(), 1);
    /// assert!(reference_sequences.contains_key(&b"sq0"[..]));
    /// # Ok::<(), Box<dyn std::error::Error>>(())
    /// ```
    pub fn reference_sequences_mut(&mut self) -> &mut ReferenceSequences {
        &mut self.reference_sequences
    }

    /// Returns the SAM header read groups.
    ///
    /// # Examples
    ///
    /// ```
    /// use pasta::sam::{
    ///     self as sam,
    ///     header::record::value::{map::ReadGroup, Map},
    /// };
    ///
    /// let header = sam::Header::builder()
    ///     .add_read_group("rg0", Map::<ReadGroup>::default())
    ///     .build();
    ///
    /// let read_groups = header.read_groups();
    /// assert_eq!(read_groups.len(), 1);
    /// assert!(read_groups.contains_key(&b"rg0"[..]));
    /// ```
    pub fn read_groups(&self) -> &ReadGroups {
        &self.read_groups
    }

    /// Returns a mutable reference to the SAM header read groups.
    ///
    /// # Examples
    ///
    /// ```
    /// use pasta::sam::{
    ///     self as sam,
    ///     header::record::value::{map::ReadGroup, Map},
    /// };
    ///
    /// let mut header = sam::Header::default();
    /// assert!(header.read_groups().is_empty());
    ///
    /// let read_group = Map::<ReadGroup>::default();
    /// header.read_groups_mut().insert(String::from("rg0").into(), read_group);
    ///
    /// let read_groups = header.read_groups();
    /// assert_eq!(read_groups.len(), 1);
    /// assert!(read_groups.contains_key(&b"rg0"[..]));
    /// ```
    pub fn read_groups_mut(&mut self) -> &mut ReadGroups {
        &mut self.read_groups
    }

    /// Returns the SAM header programs.
    ///
    /// # Examples
    ///
    /// ```
    /// use pasta::sam::{self as sam, header::record::value::{map::Program, Map}};
    ///
    /// let program = Map::<Program>::default();
    /// let header = sam::Header::builder().add_program("noodles-sam", program).build();
    ///
    /// let programs = header.programs();
    /// assert_eq!(programs.as_ref().len(), 1);
    /// assert!(programs.as_ref().contains_key(&b"noodles-sam"[..]));
    /// ```
    pub fn programs(&self) -> &Programs {
        &self.programs
    }

    /// Returns a mutable reference to the SAM header programs.
    ///
    /// # Examples
    ///
    /// ```
    /// use pasta::sam::{self as sam, header::record::value::{map::Program, Map}};
    ///
    /// let mut header = sam::Header::default();
    ///
    /// let program = Map::<Program>::default();
    /// header
    ///     .programs_mut()
    ///     .as_mut()
    ///     .insert(String::from("noodles-sam").into(), program);
    ///
    /// let programs = header.programs();
    /// assert_eq!(programs.as_ref().len(), 1);
    /// assert!(programs.as_ref().contains_key(&b"noodles-sam"[..]));
    /// ```
    pub fn programs_mut(&mut self) -> &mut Programs {
        &mut self.programs
    }

    /// Returns the SAM header comments.
    ///
    /// # Examples
    ///
    /// ```
    /// use bstr::BString;
    /// use pasta::sam;
    /// let header = sam::Header::builder().add_comment("noodles-sam").build();
    /// let comments = header.comments();
    /// assert_eq!(header.comments(), [BString::from("noodles-sam")]);
    /// ```
    pub fn comments(&self) -> &[BString] {
        &self.comments
    }

    /// Returns a mutable reference to the SAM header comments.
    ///
    /// To simply append a comment record, consider using [`Self::add_comment`] instead.
    ///
    /// # Examples
    ///
    /// ```
    /// use bstr::BString;
    /// use pasta::sam;
    /// let mut header = sam::Header::default();
    /// header.comments_mut().push(BString::from("noodles-sam"));
    /// assert_eq!(header.comments(), [BString::from("noodles-sam")]);
    /// ```
    pub fn comments_mut(&mut self) -> &mut Vec<BString> {
        &mut self.comments
    }

    /// Adds a comment.
    ///
    /// # Examples
    ///
    /// ```
    /// use bstr::BString;
    /// use pasta::sam;
    /// let mut header = sam::Header::default();
    /// header.add_comment("noodles-sam");
    /// assert_eq!(header.comments(), [BString::from("noodles-sam")]);
    /// ```
    pub fn add_comment<C>(&mut self, comment: C)
    where
        C: Into<BString>,
    {
        self.comments.push(comment.into());
    }

    /// Returns whether there are no records in this SAM header.
    ///
    /// # Examples
    ///
    /// ```
    /// use pasta::sam;
    ///
    /// let header = sam::Header::default();
    /// assert!(header.is_empty());
    ///
    /// let header = sam::Header::builder().add_comment("noodles-sam").build();
    /// assert!(!header.is_empty());
    /// ```
    pub fn is_empty(&self) -> bool {
        self.header.is_none()
            && self.reference_sequences.is_empty()
            && self.read_groups.is_empty()
            && self.programs.as_ref().is_empty()
            && self.comments.is_empty()
    }

    /// Removes all records from the header.
    ///
    /// # Examples
    ///
    /// ```
    /// use pasta::sam;
    ///
    /// let mut header = sam::Header::builder().add_comment("ndls").build();
    /// assert!(!header.is_empty());
    ///
    /// header.clear();
    /// assert!(header.is_empty());
    /// ```
    pub fn clear(&mut self) {
        self.header.take();
        self.reference_sequences.clear();
        self.read_groups.clear();
        self.programs.as_mut().clear();
        self.comments.clear();
        self.layout = None;
    }

    /// Adds an `@PG` record for a program, as htslib's `sam_hdr_add_pg` does.
    ///
    /// This is [`Programs::add`], taking into account every `@PG` line of a header that was
    /// read, including those with the same ID as a later one.
    ///
    /// # Errors
    ///
    /// This returns an error if `map` has a previous program ID (`PP`) that is not in the
    /// header.
    ///
    /// # Examples
    ///
    /// ```
    /// use pasta::sam::{
    ///     self as sam,
    ///     header::record::value::{map::{program::tag, Program}, Map},
    /// };
    ///
    /// let mut parser = sam::header::Parser::lenient();
    /// parser.parse_partial(b"@PG\tID:bwa\tPN:bwa")?;
    /// parser.parse_partial(b"@PG\tID:bwa\tPN:bwa\tVN:2")?;
    /// let mut header = parser.finish();
    ///
    /// header.add_program_record("ntools", Map::<Program>::default())?;
    ///
    /// let ntools = &header.programs().as_ref()[&b"ntools"[..]];
    /// let previous_program_id = ntools.other_fields().get(&tag::PREVIOUS_PROGRAM_ID);
    /// assert_eq!(previous_program_id.map(|id| id.as_slice()), Some(&b"bwa"[..]));
    /// # Ok::<_, Box<dyn std::error::Error>>(())
    /// ```
    pub fn add_program_record<P>(&mut self, name: P, map: Map<Program>) -> io::Result<()>
    where
        P: Into<BString>,
    {
        let lines = programs::lines(self);
        let ids: Vec<_> = programs::add(&lines, name.into(), map)?;

        for (id, map) in ids {
            self.programs.as_mut().insert(id, map);
        }

        Ok(())
    }

    /// Adds a read group at the end of the header, as htslib's `sam_hdr_add_lines` adds an
    /// `@RG` line.
    ///
    /// Inserting into [`Self::read_groups_mut`] adds a read group after the last `@RG` line of a
    /// header that was read, as htslib's `sam_hdr_add_line` does. This adds it after every other
    /// line instead, and replaces a read group with the same ID, whose line is removed.
    ///
    /// # Examples
    ///
    /// ```
    /// use pasta::sam::{self as sam, header::record::value::{map::ReadGroup, Map}};
    ///
    /// let mut parser = sam::header::Parser::lenient();
    /// parser.parse_partial(b"@RG\tID:rg0")?;
    /// parser.parse_partial(b"@PG\tID:pg0")?;
    /// let mut header = parser.finish();
    ///
    /// header.push_read_group("rg1", Map::<ReadGroup>::default())?;
    ///
    /// let mut writer = sam::io::Writer::new(Vec::new());
    /// writer.write_header(&header)?;
    /// assert_eq!(writer.into_inner(), b"@RG\tID:rg0\n@PG\tID:pg0\n@RG\tID:rg1\n");
    /// # Ok::<_, Box<dyn std::error::Error>>(())
    /// ```
    pub fn push_read_group<I>(&mut self, id: I, map: Map<ReadGroup>) -> io::Result<()>
    where
        I: Into<BString>,
    {
        let id = id.into();
        self.read_groups.shift_remove(&id);

        if let Some(layout) = self.layout.as_mut() {
            let line = Self::builder()
                .add_read_group(id.clone(), map.clone())
                .build();
            let mut writer = crate::sam::io::Writer::new(Vec::new());
            writer.write_header(&line)?;
            let src = writer.into_inner();
            let src = src.strip_suffix(b"\n").unwrap_or(&src);

            Arc::make_mut(layout).push_record(layout::Kind::ReadGroup(id.clone()), src);
        }

        self.read_groups.insert(id, map);

        Ok(())
    }

    /// Removes every read group but the one with the given ID, as htslib's
    /// `sam_hdr_remove_except` does: this also removes the later `@RG` lines with that ID, which
    /// a header that was read keeps.
    ///
    /// # Examples
    ///
    /// ```
    /// use pasta::sam;
    ///
    /// let mut parser = sam::header::Parser::lenient();
    /// parser.parse_partial(b"@RG\tID:rg0\tSM:a")?;
    /// parser.parse_partial(b"@RG\tID:rg1")?;
    /// parser.parse_partial(b"@RG\tID:rg0\tSM:b")?;
    /// let mut header = parser.finish();
    ///
    /// header.retain_read_group(b"rg0");
    ///
    /// let mut writer = sam::io::Writer::new(Vec::new());
    /// writer.write_header(&header)?;
    /// assert_eq!(writer.into_inner(), b"@RG\tID:rg0\tSM:a\n");
    /// # Ok::<_, Box<dyn std::error::Error>>(())
    /// ```
    pub fn retain_read_group(&mut self, id: &[u8]) {
        self.read_groups.retain(|rg_id, _| rg_id == id);

        if let Some(layout) = self.layout.as_mut() {
            Arc::make_mut(layout).remove_other_lines(*b"RG");
        }
    }

    pub(crate) fn layout(&self) -> Option<&Layout> {
        self.layout.as_deref()
    }

    pub(crate) fn set_layout(&mut self, layout: Option<Layout>) {
        self.layout = layout.map(Arc::new);
    }
}

impl FromStr for Header {
    type Err = ParseError;

    /// Parses a raw SAM header.
    ///
    /// # Examples
    ///
    /// ```
    /// use pasta::sam;
    ///
    /// let s = "\
    /// @HD\tVN:1.6\tSO:coordinate
    /// @SQ\tSN:sq0\tLN:8
    /// @SQ\tSN:sq1\tLN:13
    /// ";
    ///
    /// let header: sam::Header = s.parse()?;
    ///
    /// assert!(header.header().is_some());
    /// assert_eq!(header.reference_sequences().len(), 2);
    /// assert!(header.read_groups().is_empty());
    /// assert!(header.programs().as_ref().is_empty());
    /// assert!(header.comments().is_empty());
    /// # Ok::<(), sam::header::ParseError>(())
    /// ```
    fn from_str(s: &str) -> Result<Self, Self::Err> {
        parser::parse(s)
    }
}
