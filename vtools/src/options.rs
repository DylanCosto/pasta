//! Options whose meaning depends on where they are on the command line.

use clap::ArgMatches;
use pasta::variant::readlist::ListSource;

/// The command-line position of an option's last occurrence.
pub(crate) fn last_index(matches: &ArgMatches, id: &str) -> Option<usize> {
    matches.indices_of(id).and_then(Iterator::max)
}

/// The value of whichever of an inline option and its file form (`-r`/`-R`, `-s`/`-S`, ...)
/// was given last, with where to read the list from. Each is given as its argument ID and
/// value.
pub(crate) fn later<'a>(
    matches: &ArgMatches,
    (inline_id, inline): (&str, Option<&'a str>),
    (file_id, file): (&str, Option<&'a str>),
) -> Option<(&'a str, ListSource)> {
    match (inline, file) {
        (Some(x), Some(_)) if last_index(matches, inline_id) > last_index(matches, file_id) => {
            Some((x, ListSource::Inline))
        }
        (_, Some(y)) => Some((y, ListSource::File)),
        (Some(x), None) => Some((x, ListSource::Inline)),
        (None, None) => None,
    }
}

/// [`later`] for the commands whose list is a file whenever the file form was given at all,
/// even before the inline one (bcftools sets their `*_is_file` flag for good).
pub(crate) fn later_file_sticky<'a>(
    matches: &ArgMatches,
    inline: (&str, Option<&'a str>),
    file: (&str, Option<&'a str>),
) -> Option<(&'a str, ListSource)> {
    let (value, _) = later(matches, inline, file)?;
    let source = if file.1.is_some() {
        ListSource::File
    } else {
        ListSource::Inline
    };
    Some((value, source))
}
