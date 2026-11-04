use std::{
    collections::{HashMap, HashSet},
    io,
};

use bstr::{BStr, BString, ByteSlice};
use indexmap::{IndexMap, IndexSet};

use super::{
    Header,
    layout::{self, Kind},
    record::value::map::{Map, OtherFields, Program, program::tag},
};
use crate::variant::log;

const PROGRAM: [u8; 2] = *b"PG";
const ID: [u8; 2] = *b"ID";
const PREVIOUS_PROGRAM_ID: [u8; 2] = *b"PP";

type Inner = IndexMap<BString, Map<Program>>;

#[allow(clippy::tabs_in_doc_comments)]
/// SAM header programs.
///
/// SAM header programs are header records that form program chains. A program chain is a directed
/// acyclic graph (DAG) starting at a root program and ending at a leaf program. Program edges are
/// directed forward from a parent program to its linked program using the previous program ID
/// (`PP`) field in the child program.
///
/// For example, take the following program records:
///
/// ```text
/// @PG	ID:pg0
/// @PG	ID:pg1	PP:pg0
/// ```
///
/// This creates the program chain `pg0 -> pg1`.
///
/// There can exist more than one program chain in the SAM header programs, e.g.,
///
/// ```text
/// @PG	ID:pg0
/// @PG	ID:pg1
/// @PG	ID:pg2
/// @PG	ID:pg3	PP:pg0
/// @PG	ID:pg4	PP:pg1
/// @PG	ID:pg5	PP:pg4
/// ```
///
/// This creates the program chains `pg0 -> pg3`, `pg1 -> pg4 -> pg5`, and `pg2`.
#[derive(Clone, Debug, Default, Eq, PartialEq)]
pub struct Programs(Inner);

impl Programs {
    /// Adds a program, as htslib's `sam_hdr_add_pg` does.
    ///
    /// The program is added at the end of each program chain: for each program that has a
    /// previous program (`PP`) and is not itself a previous program, a record follows it. If no
    /// program has a previous program, a record follows the last program. A record gets the
    /// ID `id_prefix` or, if that is taken, the first free `<id_prefix>.<n>` (counting from 1),
    /// and a program name (`PN`) of `id_prefix` unless `map` has one.
    ///
    /// If `map` has a previous program ID, a single record links to that program instead.
    ///
    /// Note that [`crate::sam::Header::add_program_record`] also accounts for the `@PG` lines of
    /// a header that was read that are not in the programs (those with the same ID as a later
    /// one).
    ///
    /// # Errors
    ///
    /// This returns an error if `map` has a previous program ID that is not a program.
    ///
    /// # Examples
    ///
    /// ```
    /// use pasta::sam::{
    ///     self as sam,
    ///     header::record::value::{map::{program::tag, Program}, Map},
    /// };
    ///
    /// let mut header = sam::Header::default();
    /// let programs = header.programs_mut();
    ///
    /// programs.add("pg0", Map::default())?;
    /// programs.add("pg1", Map::default())?;
    ///
    /// let expected = sam::Header::builder()
    ///     .add_program("pg0", Map::builder().insert(tag::NAME, "pg0").build()?)
    ///     .add_program(
    ///         "pg1",
    ///         Map::builder()
    ///             .insert(tag::NAME, "pg1")
    ///             .insert(tag::PREVIOUS_PROGRAM_ID, "pg0")
    ///             .build()?,
    ///     )
    ///     .build();
    ///
    /// assert_eq!(programs, expected.programs());
    /// # Ok::<_, Box<dyn std::error::Error>>(())
    /// ```
    pub fn add<P>(&mut self, id_prefix: P, map: Map<Program>) -> io::Result<()>
    where
        P: Into<BString>,
    {
        let lines: Vec<_> = self
            .0
            .iter()
            .map(|(id, map)| Line::new(id.as_ref(), previous_program_id_field(map)))
            .collect();

        let records = add(&lines, id_prefix.into(), map)?;
        self.0.extend(records);

        Ok(())
    }

    /// Returns an iterator over root programs.
    ///
    /// A root program is a first program of a program chain.
    ///
    /// # Examples
    ///
    /// ```
    /// use pasta::sam::{
    ///     self as sam,
    ///     header::record::value::{map::{program::tag, Program}, Map},
    /// };
    ///
    /// let header = sam::Header::builder()
    ///     .add_program("pg0", Map::default())
    ///     .add_program("pg1", Map::builder().insert(tag::PREVIOUS_PROGRAM_ID, "pg0").build()?)
    ///     .add_program("pg2", Map::builder().insert(tag::PREVIOUS_PROGRAM_ID, "pg1").build()?)
    ///     .add_program("pg3", Map::default())
    ///     .build();
    ///
    /// let mut roots = header.programs().roots();
    /// assert_eq!(roots.next().map(|(id, _)| id.as_ref()), Some(&b"pg0"[..]));
    /// assert_eq!(roots.next().map(|(id, _)| id.as_ref()), Some(&b"pg3"[..]));
    /// assert!(roots.next().is_none());
    /// # Ok::<_, sam::header::record::value::map::builder::BuildError>(())
    /// ```
    pub fn roots(&self) -> impl Iterator<Item = (&BStr, &Map<Program>)> {
        self.0
            .iter()
            .filter(|(id, map)| previous_program_id(id.as_ref(), map).is_none())
            .map(|(id, map)| (id.as_ref(), map))
    }

    /// Returns an iterator over leaf programs.
    ///
    /// A leaf program is the last program of a program chain.
    ///
    /// # Errors
    ///
    /// This returns an `io::Error` if any program chain has a cycle. As in htslib, a program whose
    /// previous program ID is its own ID is not a cycle: the link is ignored with a warning, and the
    /// program starts a chain.
    ///
    /// # Examples
    ///
    /// ```
    /// use pasta::sam::{
    ///     self as sam,
    ///     header::record::value::{map::{program::tag, Program}, Map},
    /// };
    ///
    /// let header = sam::Header::builder()
    ///     .add_program("pg0", Map::default())
    ///     .add_program("pg1", Map::builder().insert(tag::PREVIOUS_PROGRAM_ID, "pg0").build()?)
    ///     .add_program("pg2", Map::builder().insert(tag::PREVIOUS_PROGRAM_ID, "pg1").build()?)
    ///     .add_program("pg3", Map::default())
    ///     .build();
    ///
    /// let mut leaves = header.programs().leaves()?;
    /// assert_eq!(leaves.next().map(|(id, _)| id.as_ref()), Some(&b"pg3"[..]));
    /// assert_eq!(leaves.next().map(|(id, _)| id.as_ref()), Some(&b"pg2"[..]));
    /// assert!(leaves.next().is_none());
    /// # Ok::<_, Box<dyn std::error::Error>>(())
    /// ```
    pub fn leaves(&self) -> io::Result<impl Iterator<Item = (&BStr, &Map<Program>)>> {
        let mut ids: IndexSet<_> = self.0.keys().collect();

        for (id, map) in &self.0 {
            if let Some(previous_program_id) = map.other_fields().get(&tag::PREVIOUS_PROGRAM_ID)
                && previous_program_id == id
            {
                log::warn(
                    "sam_hdr_link_pg",
                    format_args!("PG line with ID:{id} has a PP link to itself"),
                );
            } else if let Some(previous_program_id) =
                map.other_fields().get(&tag::PREVIOUS_PROGRAM_ID)
            {
                if has_cycle(&self.0, previous_program_id.as_ref(), id.as_ref()) {
                    return Err(io::Error::new(io::ErrorKind::InvalidData, "cycle detected"));
                }

                ids.swap_remove(previous_program_id);
            }
        }

        Ok(ids.into_iter().map(|id| {
            // SAFETY: `id` is guaranteed to be in the set of keys.
            self.0
                .get_key_value(id)
                .map(|(i, map)| (i.as_ref(), map))
                .unwrap()
        }))
    }
}

/// Returns a program's previous program ID, unless it is the program's own ID, a link htslib
/// ignores (`sam_hdr_link_pg`).
fn previous_program_id<'a>(id: &BStr, map: &'a Map<Program>) -> Option<&'a BStr> {
    map.other_fields()
        .get(&tag::PREVIOUS_PROGRAM_ID)
        .map(|previous_program_id| previous_program_id.as_ref())
        .filter(|&previous_program_id| previous_program_id != id)
}

/// Returns whether following the previous program IDs from `parent_id` leads back to `node_id`.
///
/// A previous program ID that is not in the graph ends the chain (HTSlib warns about such links
/// and keeps the line). The walk is bounded by the graph's size, so a cycle that does not pass
/// through `node_id` also ends it; that cycle is found from one of its own nodes.
fn has_cycle<'a>(graph: &'a Inner, mut parent_id: &'a BStr, node_id: &'a BStr) -> bool {
    for _ in 0..graph.len() {
        let Some(node) = graph.get(parent_id) else {
            return false;
        };

        if let Some(previous_program_id) = node.other_fields().get(&tag::PREVIOUS_PROGRAM_ID) {
            parent_id = previous_program_id.as_ref();

            if parent_id == node_id {
                return true;
            }
        } else {
            return false;
        }
    }

    false
}

/// An `@PG` line: its ID and previous program ID (`PP`), as htslib reads them.
#[derive(Clone, Copy, Debug)]
pub(super) struct Line<'a> {
    id: &'a BStr,
    previous_program_id: Option<&'a BStr>,
}

impl<'a> Line<'a> {
    fn new(id: &'a BStr, previous_program_id: Option<&'a BStr>) -> Self {
        Self {
            id,
            previous_program_id,
        }
    }
}

fn previous_program_id_field(map: &Map<Program>) -> Option<&BStr> {
    map.other_fields()
        .get(&tag::PREVIOUS_PROGRAM_ID)
        .map(|previous_program_id| previous_program_id.as_ref())
}

/// Returns a header's `@PG` lines in the order they are written.
///
/// For a header that was read, these are its `@PG` lines as read, including those that a later
/// line with the same ID replaces in the programs, followed by programs added since.
pub(super) fn lines(header: &Header) -> Vec<Line<'_>> {
    let programs = &header.programs().0;

    let Some(layout) = header.layout() else {
        return programs
            .iter()
            .map(|(id, map)| Line::new(id.as_ref(), previous_program_id_field(map)))
            .collect();
    };

    let mut lines = Vec::new();
    let mut in_layout = HashSet::new();

    for line in layout.lines() {
        match line.kind() {
            Kind::Program(id) => {
                in_layout.insert(id.as_bstr());

                if let Some((id, map)) = programs.get_key_value(id) {
                    lines.push(Line::new(id.as_ref(), previous_program_id_field(map)));
                }
            }
            Kind::Other => {
                let src = layout.line(line);

                if layout::record_type(src) == PROGRAM
                    && !layout::is_removed_duplicate(header, src)
                    && let Some(id) = layout::field(src, ID)
                {
                    let previous_program_id = layout::field(src, PREVIOUS_PROGRAM_ID);
                    lines.push(Line::new(
                        id.as_bstr(),
                        previous_program_id.map(|v| v.as_bstr()),
                    ));
                }
            }
            _ => {}
        }
    }

    lines.extend(
        programs
            .iter()
            .filter(|(id, _)| !in_layout.contains(id.as_bstr()))
            .map(|(id, map)| Line::new(id.as_ref(), previous_program_id_field(map))),
    );

    lines
}

/// Returns the records that adding a program to a header with the given `@PG` lines adds, as
/// htslib's `sam_hdr_add_pg` does.
pub(super) fn add(
    lines: &[Line<'_>],
    name: BString,
    map: Map<Program>,
) -> io::Result<Vec<(BString, Map<Program>)>> {
    let mut ids: HashSet<&[u8]> = lines.iter().map(|line| line.id.as_bytes()).collect();

    let previous_program_ids: Vec<Option<BString>> = if let Some(previous_program_id) =
        previous_program_id_field(&map)
    {
        if !ids.contains(previous_program_id.as_bytes()) {
            return Err(io::Error::new(
                io::ErrorKind::InvalidInput,
                format!("Header @PG ID:{previous_program_id} referred to by PP tag not present"),
            ));
        }

        vec![Some(previous_program_id.into())]
    } else {
        let ends = chain_ends(lines);

        // A cycle through every line leaves no ends.
        if ends.is_empty() {
            vec![None]
        } else {
            ends.into_iter().map(|i| Some(lines[i].id.into())).collect()
        }
    };

    // `sam_hdr_pg_id` counts on from where it left off.
    let mut n = 1;
    let mut new_ids = Vec::with_capacity(previous_program_ids.len());

    for _ in &previous_program_ids {
        let id = if ids.contains(name.as_slice()) || new_ids.contains(&name) {
            loop {
                let candidate = BString::from(format!("{name}.{n}"));
                n += 1;

                if !ids.contains(candidate.as_slice()) && !new_ids.contains(&candidate) {
                    break candidate;
                }
            }
        } else {
            name.clone()
        };

        new_ids.push(id);
    }

    ids.clear();

    let program_name = map
        .other_fields()
        .get(&tag::NAME)
        .cloned()
        .unwrap_or_else(|| name.clone());

    Ok(new_ids
        .into_iter()
        .zip(previous_program_ids)
        .map(|(id, previous_program_id)| {
            // htslib writes `ID`, `PN`, and `PP` first.
            let mut record = map.clone();
            let other_fields = record.other_fields_mut();
            other_fields.shift_remove(&tag::NAME);
            other_fields.shift_remove(&tag::PREVIOUS_PROGRAM_ID);

            let mut fields = OtherFields::with_capacity(other_fields.len() + 2);
            fields.insert(tag::NAME, program_name.clone());

            if let Some(previous_program_id) = previous_program_id {
                fields.insert(tag::PREVIOUS_PROGRAM_ID, previous_program_id);
            }

            fields.extend(other_fields.drain(..));
            *other_fields = fields;

            (id, record)
        })
        .collect())
}

/// Returns the indices of the lines that end program chains, as htslib's `sam_hdr_link_pg`
/// finds them.
///
/// A `PP` links to the last line with that ID. The chain ends are the lines that have a linked
/// `PP` and that no `PP` links to. If there are none, the end is the last line that no `PP`
/// links to. Links to a missing program or to the line itself are ignored with a warning.
fn chain_ends(lines: &[Line<'_>]) -> Vec<usize> {
    let mut indices: HashMap<&BStr, usize> = HashMap::with_capacity(lines.len());

    for (i, line) in lines.iter().enumerate() {
        indices.insert(line.id, i);
    }

    let mut is_end = vec![true; lines.len()];
    let mut chain_sizes = vec![0usize; lines.len()];

    for (i, line) in lines.iter().enumerate() {
        let Some(previous_program_id) = line.previous_program_id else {
            continue;
        };

        match indices.get(previous_program_id) {
            None => log::warn(
                "sam_hdr_link_pg",
                format_args!(
                    "PG line with ID:{} has a PP link to missing program '{previous_program_id}'",
                    line.id
                ),
            ),
            Some(&j) if j == i => log::warn(
                "sam_hdr_link_pg",
                format_args!("PG line with ID:{} has a PP link to itself", line.id),
            ),
            Some(&j) => {
                is_end[j] = false;
                chain_sizes[i] = chain_sizes[j] + 1;
            }
        }
    }

    let mut ends = Vec::new();
    let mut last_end = None;

    for (i, &end) in is_end.iter().enumerate() {
        if end {
            last_end = Some(i);

            if chain_sizes[i] > 0 {
                ends.push(i);
            }
        }
    }

    if ends.is_empty() {
        ends.extend(last_end);
    }

    ends
}

impl AsRef<Inner> for Programs {
    fn as_ref(&self) -> &Inner {
        &self.0
    }
}

impl AsMut<Inner> for Programs {
    fn as_mut(&mut self) -> &mut Inner {
        &mut self.0
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::sam::Header;

    #[test]
    fn test_add() -> Result<(), Box<dyn std::error::Error>> {
        fn program(name: &str, previous_program_id: Option<&str>) -> Map<Program> {
            let mut map = Map::<Program>::default();
            map.other_fields_mut().insert(tag::NAME, name.into());

            if let Some(id) = previous_program_id {
                map.other_fields_mut()
                    .insert(tag::PREVIOUS_PROGRAM_ID, id.into());
            }

            map
        }

        let mut programs = Programs::default();
        assert!(programs.as_ref().is_empty());

        programs.add("pg0", Map::default())?;
        let expected = Programs([("pg0".into(), program("pg0", None))].into_iter().collect());
        assert_eq!(programs, expected);

        // Without chains, a program follows the last program.
        programs.add("pg1", Map::default())?;
        assert_eq!(programs.as_ref()[&b"pg1"[..]], program("pg1", Some("pg0")));

        // Only chain ends get a program: `pg1` ends a chain, `pg2` is not in one.
        programs.as_mut().insert("pg2".into(), Map::default());
        programs.add("pg3", Map::default())?;
        assert_eq!(programs.as_ref().len(), 4);
        assert_eq!(programs.as_ref()[&b"pg3"[..]], program("pg3", Some("pg1")));

        // Two chain ends get a program each, the second with a new ID.
        let mut programs = Programs(
            [
                ("a".into(), Map::default()),
                ("b".into(), program("b", Some("a"))),
                ("c".into(), program("c", Some("a"))),
                ("d".into(), Map::default()),
            ]
            .into_iter()
            .collect(),
        );
        programs.add("b", Map::default())?;
        let ids: Vec<_> = programs.as_ref().keys().collect();
        assert_eq!(ids, ["a", "b", "c", "d", "b.1", "b.2"]);
        assert_eq!(programs.as_ref()[&b"b.1"[..]], program("b", Some("b")));
        assert_eq!(programs.as_ref()[&b"b.2"[..]], program("b", Some("c")));

        // A cycle is not an error. If every line is in one, the program starts a chain.
        let mut programs = Programs(
            [
                ("a".into(), program("a", Some("b"))),
                ("b".into(), program("b", Some("a"))),
            ]
            .into_iter()
            .collect(),
        );
        programs.add("c", Map::default())?;
        assert_eq!(programs.as_ref()[&b"c"[..]], program("c", None));

        // A given previous program ID is used, and must exist.
        programs.add("d", program("d", Some("a")))?;
        assert_eq!(programs.as_ref()[&b"d"[..]], program("d", Some("a")));
        assert!(programs.add("e", program("e", Some("nosuch"))).is_err());

        Ok(())
    }

    #[test]
    fn test_leaves() -> Result<(), Box<dyn std::error::Error>> {
        let header = Header::builder()
            .add_program("pg0", Map::default())
            .add_program("pg1", Map::default())
            .add_program("pg2", Map::default())
            .add_program(
                "pg3",
                Map::builder()
                    .insert(tag::PREVIOUS_PROGRAM_ID, "pg0")
                    .build()?,
            )
            .add_program(
                "pg4",
                Map::builder()
                    .insert(tag::PREVIOUS_PROGRAM_ID, "pg1")
                    .build()?,
            )
            .add_program(
                "pg5",
                Map::builder()
                    .insert(tag::PREVIOUS_PROGRAM_ID, "pg4")
                    .build()?,
            )
            .build();

        let mut leaves = header.programs().leaves()?;
        assert_eq!(leaves.next().map(|(id, _)| id.as_ref()), Some(&b"pg5"[..]));
        assert_eq!(leaves.next().map(|(id, _)| id.as_ref()), Some(&b"pg3"[..]));
        assert_eq!(leaves.next().map(|(id, _)| id.as_ref()), Some(&b"pg2"[..]));
        assert!(leaves.next().is_none());

        Ok(())
    }

    #[test]
    fn test_leaves_with_cycle()
    -> Result<(), crate::sam::header::record::value::map::builder::BuildError> {
        let header = Header::builder()
            .add_program(
                "pg0",
                Map::builder()
                    .insert(tag::PREVIOUS_PROGRAM_ID, "pg1")
                    .build()?,
            )
            .add_program(
                "pg1",
                Map::builder()
                    .insert(tag::PREVIOUS_PROGRAM_ID, "pg0")
                    .build()?,
            )
            .build();

        assert!(matches!(
            header.programs().leaves(),
            Err(e) if e.kind() == io::ErrorKind::InvalidData
        ));

        Ok(())
    }

    #[test]
    fn test_leaves_and_roots_with_a_self_link() -> Result<(), Box<dyn std::error::Error>> {
        // A re-run program, e.g., `@PG ID:bwa` then `@PG ID:bwa PP:bwa`, read leniently. htslib
        // looks up `PP:bwa` as the last `bwa`, i.e., the program itself, and ignores the link.
        let header = Header::builder()
            .add_program("bwa.1", Map::default())
            .add_program(
                "bwa",
                Map::builder()
                    .insert(tag::PREVIOUS_PROGRAM_ID, "bwa")
                    .build()?,
            )
            .build();

        let leaves: Vec<_> = header.programs().leaves()?.map(|(id, _)| id).collect();
        assert_eq!(leaves, [&b"bwa.1"[..], &b"bwa"[..]]);

        let roots: Vec<_> = header.programs().roots().map(|(id, _)| id).collect();
        assert_eq!(roots, [&b"bwa.1"[..], &b"bwa"[..]]);

        // htslib ignores the link, so there are no chains, and the new program follows the last.
        let mut programs = header.programs().clone();
        programs.add("pg", Map::default())?;
        assert_eq!(
            programs.as_ref()[&b"pg"[..]]
                .other_fields()
                .get(&tag::PREVIOUS_PROGRAM_ID),
            Some(&BString::from("bwa"))
        );

        Ok(())
    }

    #[test]
    fn test_leaves_with_a_missing_previous_program_or_a_downstream_cycle()
    -> Result<(), Box<dyn std::error::Error>> {
        let header = Header::builder()
            .add_program(
                "pg0",
                Map::builder()
                    .insert(tag::PREVIOUS_PROGRAM_ID, "missing")
                    .build()?,
            )
            .build();

        let leaves: Vec<_> = header.programs().leaves()?.map(|(id, _)| id).collect();
        assert_eq!(leaves, [&b"pg0"[..]]);

        // pg0 -> pg1 <-> pg2: the cycle does not pass through pg0.
        let header = Header::builder()
            .add_program(
                "pg0",
                Map::builder()
                    .insert(tag::PREVIOUS_PROGRAM_ID, "pg1")
                    .build()?,
            )
            .add_program(
                "pg1",
                Map::builder()
                    .insert(tag::PREVIOUS_PROGRAM_ID, "pg2")
                    .build()?,
            )
            .add_program(
                "pg2",
                Map::builder()
                    .insert(tag::PREVIOUS_PROGRAM_ID, "pg1")
                    .build()?,
            )
            .build();

        assert!(matches!(
            header.programs().leaves(),
            Err(e) if e.kind() == io::ErrorKind::InvalidData
        ));

        Ok(())
    }
}
