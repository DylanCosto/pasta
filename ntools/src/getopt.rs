//! GNU `getopt_long`, for commands whose options apply in command-line order or interact, so that
//! they are parsed as their C originals parse them: permuted operands, bundled short options,
//! attached or separate values, and unambiguous long option prefixes.

use anyhow::bail;

/// A short option: its character, its key, and whether it takes an argument.
pub type ShortOption<'k> = (char, &'k str, bool);

/// A long option: its name, its key, and whether it takes an argument. A key of `?` marks an
/// option the command rejects.
pub type LongOption<'k> = (&'k str, &'k str, bool);

/// The options in command-line order (as keys, with their values) and the operands.
pub type Parsed<'k> = (Vec<(&'k str, Option<String>)>, Vec<String>);

/// Parses `args` with GNU `getopt_long` rules.
///
/// # Errors
///
/// Returns an error for an unknown, ambiguous, or rejected option, a missing argument, or an
/// argument given to an option that takes none.
pub fn getopt<'k>(
    args: &[String],
    short: &[ShortOption<'k>],
    long: &[LongOption<'k>],
) -> anyhow::Result<Parsed<'k>> {
    let mut opts = Vec::new();
    let mut operands = Vec::new();
    let mut i = 0;

    while i < args.len() {
        let arg = &args[i];
        i += 1;

        if arg == "--" {
            operands.extend(args[i..].iter().cloned());
            break;
        }

        if let Some(name_value) = arg.strip_prefix("--") {
            let (name, value) = match name_value.split_once('=') {
                Some((name, value)) => (name, Some(value.to_owned())),
                None => (name_value, None),
            };

            let candidates: Vec<_> = match long.iter().find(|(n, _, _)| *n == name) {
                Some(exact) => vec![exact],
                None => long
                    .iter()
                    .filter(|(n, _, _)| n.starts_with(name))
                    .collect(),
            };

            let &(_, key, has_arg) = match candidates.as_slice() {
                [] => bail!("unrecognized option '--{name}'"),
                [first, rest @ ..]
                    if rest.iter().all(|(_, k, a)| *k == first.1 && *a == first.2) =>
                {
                    *first
                }
                _ => bail!("option '--{name}' is ambiguous"),
            };

            if key == "?" {
                bail!("unsupported option '--{name}'");
            }

            let value = if has_arg {
                match value {
                    Some(value) => Some(value),
                    None if i < args.len() => {
                        i += 1;
                        Some(args[i - 1].clone())
                    }
                    None => bail!("option '--{name}' requires an argument"),
                }
            } else if value.is_some() {
                bail!("option '--{name}' doesn't allow an argument");
            } else {
                None
            };

            opts.push((key, value));
        } else if let Some(cluster) = arg.strip_prefix('-').filter(|s| !s.is_empty()) {
            for (j, c) in cluster.char_indices() {
                let Some(&(_, key, has_arg)) = short.iter().find(|(s, _, _)| *s == c) else {
                    bail!("invalid option -- '{c}'");
                };

                if has_arg {
                    let rest = &cluster[j + c.len_utf8()..];

                    let value = if !rest.is_empty() {
                        rest.to_owned()
                    } else if i < args.len() {
                        i += 1;
                        args[i - 1].clone()
                    } else {
                        bail!("option requires an argument -- '{c}'");
                    };

                    opts.push((key, Some(value)));
                    break;
                }

                opts.push((key, None));
            }
        } else {
            operands.push(arg.clone());
        }
    }

    Ok((opts, operands))
}

/// C's `atol`: leading whitespace, an optional sign, and digits; 0 if there are none.
pub fn atol(s: &str) -> i64 {
    let s = s.trim_start_matches(|c: char| c.is_ascii_whitespace());
    let (neg, digits) = match s.as_bytes().first() {
        Some(b'-') => (true, &s[1..]),
        Some(b'+') => (false, &s[1..]),
        _ => (false, s),
    };

    let n = digits
        .bytes()
        .take_while(u8::is_ascii_digit)
        .fold(0i64, |n, d| {
            n.saturating_mul(10).saturating_add(i64::from(d - b'0'))
        });

    if neg { -n } else { n }
}

/// Parses an integer option that samtools reads with `atoi` and where a negative value means
/// the same as 0 (`view -q -5` keeps every read, `-@ -1` adds no threads).
pub fn negative_as_0<T>(s: &str) -> Result<T, String>
where
    T: std::str::FromStr + Default,
    T::Err: std::fmt::Display,
{
    match s.parse::<i64>() {
        Ok(n) if n < 0 => Ok(T::default()),
        _ => s.parse().map_err(|e: T::Err| e.to_string()),
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    const SHORT: &[ShortOption<'static>] = &[('a', "a", false), ('f', "format", true)];
    const LONG: &[LongOption<'static>] = &[
        ("format", "format", true),
        ("show-del", "show-del", true),
        ("show-ins", "show-ins", true),
        ("gone", "?", false),
    ];

    fn args(s: &str) -> Vec<String> {
        s.split_whitespace().map(String::from).collect()
    }

    #[test]
    fn test_getopt() -> anyhow::Result<()> {
        let (opts, operands) = getopt(
            &args("-aa in.bam -fpileup --show-d yes - --format=x -- -a"),
            SHORT,
            LONG,
        )?;

        assert_eq!(operands, ["in.bam", "-", "-a"]);
        assert_eq!(
            opts,
            [
                ("a", None),
                ("a", None),
                ("format", Some("pileup".into())),
                ("show-del", Some("yes".into())),
                ("format", Some("x".into())),
            ]
        );

        assert!(getopt(&args("--show 3"), SHORT, LONG).is_err());
        assert!(getopt(&args("--gone"), SHORT, LONG).is_err());
        assert!(getopt(&args("-h"), SHORT, LONG).is_err());
        assert!(getopt(&args("-f"), SHORT, LONG).is_err());
        assert!(getopt(&args("--format"), SHORT, LONG).is_err());

        Ok(())
    }

    #[test]
    fn test_atol() {
        assert_eq!(atol("120"), 120);
        assert_eq!(atol(" -5x"), -5);
        assert_eq!(atol("+7"), 7);
        assert_eq!(atol("abc"), 0);
        assert_eq!(atol(""), 0);
    }
}
