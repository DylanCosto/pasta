//! POSIX basic regular expressions (`regcomp` without `REG_EXTENDED`, as bcftools compiles
//! `~` patterns), translated to the `regex` crate's syntax. Patterns with back-references
//! (`\1` to `\9`), which no finite automaton matches, run on a memoizing backtracking
//! matcher compiled from the translation's syntax tree.

use std::{cell::RefCell, collections::HashSet, fmt::Write as _};

use regex::bytes::{Regex, RegexBuilder};
use regex_syntax::hir::{Class, Hir, HirKind, Look};

/// Whether a pattern matches letters of either case (`REG_ICASE`).
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub(crate) enum Case {
    Sensitive,
    Insensitive,
}

/// A compiled basic regular expression.
#[derive(Clone, Debug)]
pub(crate) struct Bre(Matcher);

#[derive(Clone, Debug)]
enum Matcher {
    Automaton(Regex),
    Backtrack(Box<Backtrack>),
}

impl Bre {
    /// Whether the expression matches anywhere in `text` (`regexec`).
    pub(crate) fn is_match(&self, text: &[u8]) -> bool {
        match &self.0 {
            Matcher::Automaton(re) => re.is_match(text),
            Matcher::Backtrack(b) => b.is_match(text),
        }
    }
}

/// Compiles a POSIX basic regular expression (with GNU's `\+`, `\?`, and `\|`) as glibc does
/// in the C locale: byte-wise, `.` matching any byte.
pub(crate) fn compile(pattern: &[u8], case: Case) -> Result<Bre, String> {
    let translated = translate(pattern, Backref::AnyText)?;
    let automaton = RegexBuilder::new(&translated.pattern)
        .unicode(false)
        .dot_matches_new_line(true)
        .case_insensitive(case == Case::Insensitive)
        .build()
        .map_err(|e| e.to_string())?;
    if !translated.backrefs {
        return Ok(Bre(Matcher::Automaton(automaton)));
    }
    let marked = translate(pattern, Backref::Marker)?;
    let hir = regex_syntax::ParserBuilder::new()
        .unicode(false)
        .utf8(false)
        .dot_matches_new_line(true)
        .case_insensitive(case == Case::Insensitive)
        .build()
        .parse(&marked.pattern)
        .map_err(|e| e.to_string())?;
    let node = Node::from_hir(&hir)?;
    Ok(Bre(Matcher::Backtrack(Box::new(Backtrack {
        prefilter: automaton,
        program: Program::new(&node)?,
        case,
    }))))
}

/// How [`translate`] writes a back-reference.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
enum Backref {
    /// As `(?:.*)`, which matches whatever the reference can: the translation then matches a
    /// superset of the pattern's matches.
    AnyText,
    /// As a group named [`BACKREF_PREFIX`], a running count, `_`, and the group number, for
    /// the backtracking matcher. The group holds `.*` so that the parser does not take a
    /// repetition of it for a repetition of the empty string, which it reduces to one.
    Marker,
}

/// The start of the names of the groups that mark back-references.
const BACKREF_PREFIX: &str = "backref";

/// The start of the names of the pattern's own groups in a [`Backref::Marker`] translation,
/// followed by the group's number. The parser drops a group repeated `\{0\}` times, so the
/// matcher takes each group's number from its name rather than from counting the groups.
const GROUP_PREFIX: &str = "group";

/// A pattern in the `regex` syntax.
struct Translated {
    pattern: String,
    /// Whether the pattern has back-references.
    backrefs: bool,
}

/// Appends a literal byte outside a bracket expression.
fn push_literal(out: &mut String, c: u8) {
    if c.is_ascii() && !c.is_ascii_control() {
        out.push_str(&regex::escape(&char::from(c).to_string()));
    } else {
        let _ = write!(out, "\\x{c:02X}");
    }
}

/// Appends a literal byte inside a bracket expression.
fn push_class_literal(out: &mut String, c: u8) {
    if matches!(c, b'\\' | b'[' | b']' | b'^' | b'-' | b'&' | b'~') {
        out.push('\\');
        out.push(char::from(c));
    } else if c.is_ascii() && !c.is_ascii_control() {
        out.push(char::from(c));
    } else {
        let _ = write!(out, "\\x{c:02X}");
    }
}

/// Rewrites a BRE in the `regex` syntax, refusing back-references to a group that is not
/// closed before them in the same alternative, as glibc does.
fn translate(p: &[u8], backref: Backref) -> Result<Translated, String> {
    let mut out = String::with_capacity(p.len() * 2);
    let mut i = 0;
    let mut nbackrefs = 0;
    // The groups opened and not closed (their numbers), and those closed (a bit per number,
    // glibc's `completed_bkref_map`, as it stood at the start of each open group's current
    // alternative, and as accumulated over its earlier alternatives).
    let mut ngroups: u32 = 0;
    let mut open: Vec<Alternatives> = vec![Alternatives::default()];
    let mut closed: u32 = 0;
    // Where glibc's parser expects an expression, so that a repetition operator (`*`, `\+`,
    // `\?`) is a literal character and `\{` is refused: the start of the pattern or of a group
    // or alternative, and after an anchor, which nothing repeats.
    let mut at_start = true;
    // Where `^` is an anchor: the start of the pattern or of a group or alternative.
    let mut caret_anchors = true;
    // Where the last expression starts in `out`, and whether it ends in a repetition operator.
    let mut expr_start = 0;
    let mut repeated = false;
    while i < p.len() {
        let c = p[i];
        let mut next_start = false;
        let mut next_caret_anchors = false;
        let mut next_expr_start = out.len();
        let mut next_repeated = false;
        match c {
            b'\\' => {
                let Some(&n) = p.get(i + 1) else {
                    return Err("Trailing backslash".into());
                };
                i += 1;
                match n {
                    b'(' => {
                        next_start = true;
                        next_caret_anchors = true;
                        ngroups += 1;
                        match backref {
                            Backref::AnyText => out.push('('),
                            Backref::Marker => {
                                let _ = write!(out, "(?P<{GROUP_PREFIX}{ngroups}>");
                            }
                        }
                        open.push(Alternatives {
                            group: ngroups,
                            start: next_expr_start,
                            initial: closed,
                            accumulated: 0,
                        });
                    }
                    b')' => {
                        out.push(')');
                        // An unmatched `\)` is left for the `regex` parser to refuse.
                        if open.len() > 1
                            && let Some(alt) = open.pop()
                        {
                            closed |= alt.accumulated | 1u32.checked_shl(alt.group).unwrap_or(0);
                            next_expr_start = alt.start;
                        }
                    }
                    b'|' => {
                        out.push('|');
                        next_start = true;
                        next_caret_anchors = true;
                        if let Some(alt) = open.last_mut() {
                            alt.accumulated |= closed;
                            closed = alt.initial;
                        }
                    }
                    b'{' => {
                        let close = p[i + 1..]
                            .windows(2)
                            .position(|w| w == b"\\}")
                            .ok_or("Unmatched \\{")?;
                        let interval = interval(&p[i + 1..i + 1 + close])?;
                        if at_start {
                            return Err("Invalid preceding regular expression".into());
                        }
                        push_repetition(&mut out, expr_start, repeated, &interval)?;
                        next_repeated = true;
                        i += close + 2;
                    }
                    b'+' | b'?' if !at_start => {
                        push_repetition(
                            &mut out,
                            expr_start,
                            repeated,
                            &char::from(n).to_string(),
                        )?;
                        next_repeated = true;
                    }
                    b'1'..=b'9' => {
                        let group = n - b'0';
                        if closed & (1 << group) == 0 {
                            return Err("Invalid back reference".into());
                        }
                        match backref {
                            Backref::AnyText => out.push_str("(?:.*)"),
                            Backref::Marker => {
                                let _ = write!(out, "(?P<{BACKREF_PREFIX}{nbackrefs}_{group}>.*)");
                            }
                        }
                        nbackrefs += 1;
                    }
                    b'w' | b'W' | b's' | b'S' => {
                        out.push('\\');
                        out.push(char::from(n));
                    }
                    b'b' | b'B' | b'<' | b'>' => {
                        out.push('\\');
                        out.push(char::from(n));
                        next_start = true;
                    }
                    b'`' => {
                        out.push_str("\\A");
                        next_start = true;
                    }
                    b'\'' => {
                        out.push_str("\\z");
                        next_start = true;
                    }
                    _ => push_literal(&mut out, n),
                }
                i += 1;
            }
            b'[' => {
                i = translate_bracket(p, i, &mut out)?;
            }
            b'*' if at_start => {
                out.push_str("\\*");
                i += 1;
            }
            b'*' => {
                push_repetition(&mut out, expr_start, repeated, "*")?;
                next_repeated = true;
                i += 1;
            }
            b'^' if caret_anchors => {
                out.push('^');
                next_start = true;
                i += 1;
            }
            b'$' if i + 1 == p.len()
                || p[i + 1..].starts_with(b"\\)")
                || p[i + 1..].starts_with(b"\\|") =>
            {
                out.push('$');
                next_start = true;
                i += 1;
            }
            b'.' => {
                out.push('.');
                i += 1;
            }
            _ => {
                push_literal(&mut out, c);
                i += 1;
            }
        }
        at_start = next_start;
        caret_anchors = next_caret_anchors;
        if !next_repeated {
            expr_start = next_expr_start;
        }
        repeated = next_repeated;
    }
    Ok(Translated {
        pattern: out,
        backrefs: nbackrefs > 0,
    })
}

/// glibc's `RE_DUP_MAX`, the largest count of an interval.
const DUP_MAX: u32 = 0x7fff;

/// Translates the body of an interval (`m,n` of `\{m,n\}`) as glibc's `parse_dup_op` reads it:
/// a missing lower bound is 0 (`\{,n\}` is `\{0,n\}`), and a count over [`DUP_MAX`] is refused.
fn interval(body: &[u8]) -> Result<String, String> {
    let invalid = || "Invalid content of \\{\\}".to_string();
    // glibc's tokens: `\,` and `\0` are the characters `,` and `0`; any other escape (a
    // back-reference, an operator) is not a digit.
    let mut unescaped = Vec::with_capacity(body.len());
    let mut bytes = body.iter();
    while let Some(&b) = bytes.next() {
        unescaped.push(match (b, bytes.as_slice().first()) {
            (b'\\', Some(&c @ (b',' | b'0'))) => {
                bytes.next();
                c
            }
            (b'\\', _) => {
                bytes.next();
                b'\\'
            }
            _ => b,
        });
    }
    let body = unescaped.as_slice();
    // glibc's fetch_number, which stops counting past DUP_MAX.
    let count = |digits: &[u8]| -> Result<Option<u32>, String> {
        if !digits.iter().all(u8::is_ascii_digit) {
            return Err(invalid());
        }
        Ok((!digits.is_empty()).then(|| {
            digits.iter().fold(0, |n: u32, &d| {
                (n * 10 + u32::from(d - b'0')).min(DUP_MAX + 1)
            })
        }))
    };
    let (min, max) = match body.iter().position(|&b| b == b',') {
        None => {
            let n = count(body)?.ok_or_else(invalid)?;
            (n, Some(n))
        }
        Some(comma) => (
            count(&body[..comma])?.unwrap_or(0),
            count(&body[comma + 1..])?,
        ),
    };
    if max.is_some_and(|max| min > max) {
        return Err(invalid());
    }
    if max.unwrap_or(min) > DUP_MAX {
        return Err("Regular expression too big".into());
    }
    Ok(match max {
        Some(max) if max == min => format!("{{{min}}}"),
        Some(max) => format!("{{{min},{max}}}"),
        None => format!("{{{min},}}"),
    })
}

/// Appends the repetition operator `op` to the expression at `out[start..]`. glibc applies a
/// `\+` or `\?` after a repetition to the repeated expression (`a\{2\}\?` is `(?:a{2})?`,
/// not the lazy `a{2}?`) and refuses a `*` or `\{` there.
fn push_repetition(out: &mut String, start: usize, repeated: bool, op: &str) -> Result<(), String> {
    if repeated {
        if op.starts_with(['*', '{']) {
            return Err("Invalid preceding regular expression".into());
        }
        out.insert_str(start, "(?:");
        out.push(')');
    }
    out.push_str(op);
    Ok(())
}

/// A group being parsed, for [`translate`]'s check of back-references.
#[derive(Clone, Copy, Debug, Default)]
struct Alternatives {
    /// Its number (0 for the whole pattern).
    group: u32,
    /// Where it starts in the translation.
    start: usize,
    /// The groups closed before its current alternative.
    initial: u32,
    /// The groups closed in its earlier alternatives.
    accumulated: u32,
}

/// Translates the bracket expression starting at `p[start]` (`[`), returning the index after
/// it. A backslash is literal inside.
fn translate_bracket(p: &[u8], start: usize, out: &mut String) -> Result<usize, String> {
    let mut i = start + 1;
    out.push('[');
    if p.get(i) == Some(&b'^') {
        out.push('^');
        i += 1;
    }
    // A `]` first in the list is literal.
    if p.get(i) == Some(&b']') {
        out.push_str("\\]");
        i += 1;
    }
    loop {
        let Some(&c) = p.get(i) else {
            return Err("Unmatched [, [^, [:, [., or [=".into());
        };
        match c {
            b']' => {
                out.push(']');
                return Ok(i + 1);
            }
            b'[' if matches!(p.get(i + 1), Some(b':' | b'=' | b'.')) => {
                let kind = p[i + 1];
                let close = [kind, b']'];
                let body_start = i + 2;
                let len = p[body_start..]
                    .windows(2)
                    .position(|w| w == close)
                    .ok_or("Unmatched [, [^, [:, [., or [=")?;
                let body = &p[body_start..body_start + len];
                if kind == b':' {
                    out.push_str("[:");
                    out.push_str(std::str::from_utf8(body).map_err(|e| e.to_string())?);
                    out.push_str(":]");
                } else {
                    for &b in body {
                        push_class_literal(out, b);
                    }
                }
                i = body_start + len + 2;
            }
            b'-' => {
                out.push('-');
                i += 1;
            }
            _ => {
                push_class_literal(out, c);
                i += 1;
            }
        }
    }
}

/// A pattern with back-references, matched by backtracking.
#[derive(Clone, Debug)]
struct Backtrack {
    /// The pattern with each back-reference as any text, which matches whenever the pattern
    /// does, to skip most texts that cannot match.
    prefilter: Regex,
    program: Program,
    case: Case,
}

/// The syntax tree the backtracking matcher walks.
#[derive(Clone, Debug)]
enum Node {
    Empty,
    Literal(Box<[u8]>),
    /// One byte of a set.
    Class(Box<[bool; 256]>),
    Look(Look),
    /// A group, by its index (the group number less one).
    Group(usize, Box<Node>),
    /// A back-reference to a group, by its index.
    Backref(usize),
    Concat(Vec<Node>),
    Alternation(Vec<Node>),
    Repetition {
        sub: Box<Node>,
        min: u32,
        max: Option<u32>,
    },
}

impl Node {
    /// Converts a translation with [`Backref::Marker`]s, numbering its other groups in the
    /// order they open from `*groups` on.
    fn from_hir(hir: &Hir) -> Result<Self, String> {
        let all = |subs: &[Hir]| -> Result<Vec<Node>, String> {
            subs.iter().map(Self::from_hir).collect()
        };
        Ok(match hir.kind() {
            HirKind::Empty => Self::Empty,
            HirKind::Literal(lit) => Self::Literal(lit.0.clone()),
            HirKind::Class(class) => {
                let mut set = Box::new([false; 256]);
                match class {
                    Class::Bytes(c) => {
                        for r in c.ranges() {
                            set[usize::from(r.start())..=usize::from(r.end())].fill(true);
                        }
                    }
                    Class::Unicode(c) => {
                        for r in c.ranges() {
                            for (b, member) in set.iter_mut().enumerate() {
                                let ch = u32::try_from(b).unwrap_or(u32::MAX);
                                *member |=
                                    (u32::from(r.start())..=u32::from(r.end())).contains(&ch);
                            }
                        }
                    }
                }
                Self::Class(set)
            }
            HirKind::Look(look) => match look {
                Look::Start
                | Look::End
                | Look::WordAscii
                | Look::WordAsciiNegate
                | Look::WordStartAscii
                | Look::WordEndAscii => Self::Look(*look),
                _ => return Err(format!("unsupported assertion {look:?}")),
            },
            HirKind::Repetition(rep) => Self::Repetition {
                sub: Box::new(Self::from_hir(&rep.sub)?),
                min: rep.min,
                max: rep.max,
            },
            HirKind::Capture(cap) => {
                let name = cap.name.as_deref().unwrap_or_default();
                let number = |n: &str| n.parse::<usize>().map_err(|_| format!("bad group {name}"));
                if let Some(marker) = name.strip_prefix(BACKREF_PREFIX) {
                    let group = marker.rsplit('_').next().unwrap_or_default();
                    Self::Backref(number(group)? - 1)
                } else if let Some(group) = name.strip_prefix(GROUP_PREFIX) {
                    Self::Group(number(group)? - 1, Box::new(Self::from_hir(&cap.sub)?))
                } else {
                    return Err(format!("unnamed group {name}"));
                }
            }
            HirKind::Concat(subs) => Self::Concat(all(subs)?),
            HirKind::Alternation(subs) => Self::Alternation(all(subs)?),
        })
    }
}

/// The most groups back-references can name (`\1` to `\9`).
const MAX_REFERENCED: usize = 9;

/// The most instructions a [`Program`] may have, which bounds what counted repetitions
/// expand to.
const MAX_PROGRAM: usize = 1 << 20;

/// The most bits [`Memo`] keeps a bit per instruction and position in (16 MiB); a longer
/// program or text keys every state in the hash set.
const MAX_PLAIN_BITS: usize = 1 << 27;

/// How many keyed states [`Memo`] keeps from one start position to the next, so that a
/// search from each start of a long text needs no more memory than the largest one; also
/// the most states, frames, or words of bits a thread keeps allocated between searches.
const MAX_KEPT_STATES: usize = 1 << 16;

/// A slot of a group that has not matched.
const UNSET: usize = usize::MAX;

/// One instruction of a [`Program`]. Slot `2 * g` holds where the group of index `g` last
/// started and slot `2 * g + 1` where it last ended; only groups that back-references name
/// have them.
#[derive(Clone, Debug)]
enum Inst {
    Byte(u8),
    Class(Box<[bool; 256]>),
    Look(Look),
    /// Tries the first instruction, then the second.
    Split(usize, usize),
    Jump(usize),
    Save(usize),
    /// A back-reference to a group, by its index.
    Backref(usize),
    Match,
}

/// A pattern compiled to instructions for the backtracking matcher.
#[derive(Clone, Debug)]
struct Program {
    insts: Vec<Inst>,
    /// For each instruction, the slots (a bit each) that a back-reference may read before
    /// they are saved again, the only ones on which what can match from there depends.
    live: Vec<u32>,
}

impl Program {
    fn new(node: &Node) -> Result<Self, String> {
        let mut referenced = 0;
        node.referenced(&mut referenced);
        let mut insts = Vec::new();
        emit(node, referenced, &mut insts)?;
        insts.push(Inst::Match);
        let live = liveness(&insts);
        Ok(Self { insts, live })
    }
}

impl Node {
    /// Adds a bit for each group a back-reference in the node names.
    fn referenced(&self, groups: &mut u32) {
        match self {
            Self::Empty | Self::Literal(_) | Self::Class(_) | Self::Look(_) => {}
            Self::Backref(index) => *groups |= 1 << index,
            Self::Group(_, sub) | Self::Repetition { sub, .. } => sub.referenced(groups),
            Self::Concat(subs) | Self::Alternation(subs) => {
                for sub in subs {
                    sub.referenced(groups);
                }
            }
        }
    }
}

/// Appends the instructions of `node`, saving the slots of the `referenced` groups (a bit
/// per index). Repetitions try as many iterations as possible first.
fn emit(node: &Node, referenced: u32, insts: &mut Vec<Inst>) -> Result<(), String> {
    if insts.len() > MAX_PROGRAM {
        return Err("Regular expression too big".into());
    }
    match node {
        Node::Empty => {}
        Node::Literal(lit) => insts.extend(lit.iter().map(|&b| Inst::Byte(b))),
        Node::Class(set) => insts.push(Inst::Class(set.clone())),
        Node::Look(look) => insts.push(Inst::Look(*look)),
        Node::Group(index, sub) => {
            let saved = *index < MAX_REFERENCED && referenced & (1 << index) != 0;
            if saved {
                insts.push(Inst::Save(2 * index));
            }
            emit(sub, referenced, insts)?;
            if saved {
                insts.push(Inst::Save(2 * index + 1));
            }
        }
        Node::Backref(index) => insts.push(Inst::Backref(*index)),
        Node::Concat(subs) => {
            for sub in subs {
                emit(sub, referenced, insts)?;
            }
        }
        Node::Alternation(subs) => {
            let mut jumps = Vec::new();
            let (last, rest) = subs.split_last().ok_or("empty alternation")?;
            for sub in rest {
                let split = insts.len();
                insts.push(Inst::Split(split + 1, UNSET));
                emit(sub, referenced, insts)?;
                jumps.push(insts.len());
                insts.push(Inst::Jump(UNSET));
                insts[split] = Inst::Split(split + 1, insts.len());
            }
            emit(last, referenced, insts)?;
            let end = insts.len();
            for jump in jumps {
                insts[jump] = Inst::Jump(end);
            }
        }
        Node::Repetition { sub, min, max } => {
            for _ in 0..*min {
                emit(sub, referenced, insts)?;
            }
            if let Some(max) = max {
                let mut splits = Vec::new();
                for _ in *min..*max {
                    splits.push(insts.len());
                    insts.push(Inst::Split(insts.len() + 1, UNSET));
                    emit(sub, referenced, insts)?;
                }
                let end = insts.len();
                for split in splits {
                    insts[split] = Inst::Split(split + 1, end);
                }
            } else {
                let split = insts.len();
                insts.push(Inst::Split(split + 1, UNSET));
                emit(sub, referenced, insts)?;
                insts.push(Inst::Jump(split));
                insts[split] = Inst::Split(split + 1, insts.len());
            }
        }
    }
    Ok(())
}

/// The live slots of each instruction, by the usual backward data flow: a back-reference
/// reads its group's two slots and a save writes one.
fn liveness(insts: &[Inst]) -> Vec<u32> {
    let mut live = vec![0u32; insts.len()];
    let mut changed = true;
    while changed {
        changed = false;
        for pc in (0..insts.len()).rev() {
            let at = |to: usize| live.get(to).copied().unwrap_or(0);
            let new = match insts[pc] {
                Inst::Byte(_) | Inst::Class(_) | Inst::Look(_) => at(pc + 1),
                Inst::Split(first, second) => at(first) | at(second),
                Inst::Jump(to) => at(to),
                Inst::Save(slot) => at(pc + 1) & !(1 << slot),
                Inst::Backref(index) => at(pc + 1) | 0b11 << (2 * index),
                Inst::Match => 0,
            };
            if new != live[pc] {
                live[pc] = new;
                changed = true;
            }
        }
    }
    live
}

/// The states of a search already tried at a `Split`, from which nothing matches when they
/// come again: an instruction, a position, and the values of the instruction's live slots,
/// on which alone the rest of the match depends.
#[derive(Debug, Default)]
struct Memo {
    /// A bit per instruction and position, for the states with no live slot.
    plain: Vec<u64>,
    /// The number of positions (the text's length and one).
    width: usize,
    /// The bits that hold an instruction and a position or slot value in a packed key.
    pc_bits: u32,
    pos_bits: u32,
    /// The other states whose instruction, position, and live slots' values fit in 128 bits.
    packed: HashSet<u128>,
    /// The rest, as the instruction, the position, and the live slots' values.
    keyed: HashSet<Box<[usize]>>,
    key: Vec<usize>,
}

impl Memo {
    /// Forgets every state, for a search of `program` in `text`.
    fn reset(&mut self, program: &Program, text: &[u8]) {
        self.width = text.len() + 1;
        let bits = program
            .insts
            .len()
            .checked_mul(self.width)
            .filter(|&bits| bits <= MAX_PLAIN_BITS);
        self.plain.clear();
        self.plain
            .resize(bits.map_or(0, |bits| bits.div_ceil(64)), 0);
        self.pc_bits = usize::BITS - program.insts.len().leading_zeros();
        // Positions run to the text's length; `width` stands for an unset slot.
        self.pos_bits = usize::BITS - self.width.leading_zeros();
        self.forget_keyed();
    }

    /// The number of keyed states.
    fn len(&self) -> usize {
        self.packed.len() + self.keyed.len()
    }

    /// Forgets the keyed states.
    fn forget_keyed(&mut self) {
        if !self.packed.is_empty() {
            self.packed.clear();
        }
        if !self.keyed.is_empty() {
            self.keyed.clear();
        }
    }

    /// Frees what a search of a long text left allocated.
    fn shrink(&mut self) {
        if self.plain.capacity() > MAX_KEPT_STATES {
            self.plain = Vec::new();
        }
        if self.packed.capacity() > MAX_KEPT_STATES {
            self.packed = HashSet::new();
        }
        if self.keyed.capacity() > MAX_KEPT_STATES {
            self.keyed = HashSet::new();
        }
    }

    /// Records a state, returning whether it is new.
    fn insert(&mut self, pc: usize, pos: usize, live: u32, slots: &[usize]) -> bool {
        if live == 0 && !self.plain.is_empty() {
            let bit = pc * self.width + pos;
            let word = &mut self.plain[bit / 64];
            let mask = 1 << (bit % 64);
            let new = *word & mask == 0;
            *word |= mask;
            return new;
        }
        let values = slots
            .iter()
            .enumerate()
            .filter(|&(slot, _)| live & (1 << slot) != 0)
            .map(|(_, &value)| value);
        // The instruction decides which slots are live, so packed keys do not collide.
        if self.pc_bits + self.pos_bits * (1 + live.count_ones()) <= u128::BITS {
            let mut key = pc as u128;
            let mut shift = self.pc_bits;
            for value in std::iter::once(pos).chain(values) {
                key |= (value.min(self.width) as u128) << shift;
                shift += self.pos_bits;
            }
            return self.packed.insert(key);
        }
        self.key.clear();
        self.key.extend([pc, pos]);
        self.key.extend(values);
        if self.keyed.contains(self.key.as_slice()) {
            return false;
        }
        self.keyed.insert(self.key.as_slice().into());
        true
    }
}

/// Where the matcher goes back to when a path fails.
#[derive(Clone, Copy, Debug)]
enum Frame {
    /// Tries an instruction at a position.
    Explore(usize, usize),
    /// Puts back a slot's earlier value.
    Restore(usize, usize),
}

/// What searches allocate, kept from one to the next on each thread.
#[derive(Debug, Default)]
struct Scratch {
    memo: Memo,
    stack: Vec<Frame>,
}

thread_local! {
    static SCRATCH: RefCell<Scratch> = RefCell::default();
}

/// The state of one search.
struct Search<'a> {
    program: &'a Program,
    text: &'a [u8],
    case: Case,
    slots: [usize; 2 * MAX_REFERENCED],
    stack: &'a mut Vec<Frame>,
    memo: &'a mut Memo,
}

impl Backtrack {
    fn is_match(&self, text: &[u8]) -> bool {
        if !self.prefilter.is_match(text) {
            return false;
        }
        SCRATCH.with_borrow_mut(|scratch| {
            scratch.memo.reset(&self.program, text);
            let mut search = Search {
                program: &self.program,
                text,
                case: self.case,
                slots: [UNSET; 2 * MAX_REFERENCED],
                stack: &mut scratch.stack,
                memo: &mut scratch.memo,
            };
            let found = (0..=text.len()).any(|start| search.run(start));
            scratch.memo.shrink();
            if scratch.stack.capacity() > MAX_KEPT_STATES {
                scratch.stack = Vec::new();
            }
            found
        })
    }
}

impl Search<'_> {
    /// Whether the program matches from `start`, trying every way it can, depth first with
    /// the stack on the heap. A back-reference matches the text its group last matched
    /// (glibc's behavior also for groups inside repetitions), and never a group that has not
    /// matched. Each state is tried at most once over all the starts, which bounds the time
    /// by a polynomial in the text's length.
    fn run(&mut self, start: usize) -> bool {
        if self.memo.len() > MAX_KEPT_STATES {
            self.memo.forget_keyed();
        }
        self.slots.fill(UNSET);
        self.stack.clear();
        self.stack.push(Frame::Explore(0, start));
        let program = self.program;
        while let Some(frame) = self.stack.pop() {
            let (mut pc, mut pos) = match frame {
                Frame::Restore(slot, value) => {
                    self.slots[slot] = value;
                    continue;
                }
                Frame::Explore(pc, pos) => (pc, pos),
            };
            loop {
                match &program.insts[pc] {
                    Inst::Match => return true,
                    Inst::Byte(b) => {
                        if self.text.get(pos) != Some(b) {
                            break;
                        }
                        pc += 1;
                        pos += 1;
                    }
                    Inst::Class(set) => {
                        if !self.text.get(pos).is_some_and(|&b| set[usize::from(b)]) {
                            break;
                        }
                        pc += 1;
                        pos += 1;
                    }
                    Inst::Look(look) => {
                        if !self.look(*look, pos) {
                            break;
                        }
                        pc += 1;
                    }
                    // Every loop of the program goes through a split, so remembering the
                    // states there is enough to try none twice.
                    Inst::Split(first, second) => {
                        if !self.memo.insert(pc, pos, program.live[pc], &self.slots) {
                            break;
                        }
                        self.stack.push(Frame::Explore(*second, pos));
                        pc = *first;
                    }
                    Inst::Jump(to) => pc = *to,
                    Inst::Save(slot) => {
                        self.stack.push(Frame::Restore(*slot, self.slots[*slot]));
                        self.slots[*slot] = pos;
                        pc += 1;
                    }
                    Inst::Backref(index) => {
                        let Some(end) = self.backref(*index, pos) else {
                            break;
                        };
                        pc += 1;
                        pos = end;
                    }
                }
            }
        }
        false
    }

    /// Where a back-reference to the group of `index` at `pos` ends, if it matches.
    fn backref(&self, index: usize, pos: usize) -> Option<usize> {
        let (start, end) = (self.slots[2 * index], self.slots[2 * index + 1]);
        if start == UNSET || end == UNSET {
            return None;
        }
        let there = self.text.get(start..end)?;
        let here = self.text.get(pos..pos + there.len())?;
        let same = match self.case {
            Case::Sensitive => here == there,
            Case::Insensitive => here.eq_ignore_ascii_case(there),
        };
        same.then_some(pos + there.len())
    }

    fn look(&self, look: Look, pos: usize) -> bool {
        let word = |i: Option<usize>| {
            i.and_then(|i| self.text.get(i))
                .is_some_and(|&b| b.is_ascii_alphanumeric() || b == b'_')
        };
        let (before, after) = (word(pos.checked_sub(1)), word(Some(pos)));
        match look {
            Look::Start => pos == 0,
            Look::End => pos == self.text.len(),
            Look::WordAscii => before != after,
            Look::WordAsciiNegate => before == after,
            Look::WordStartAscii => !before && after,
            Look::WordEndAscii => before && !after,
            _ => false,
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn is_match(p: &str, s: &str) -> bool {
        compile(p.as_bytes(), Case::Sensitive)
            .expect("the pattern compiles")
            .is_match(s.as_bytes())
    }

    #[test]
    fn basic_syntax() {
        assert!(is_match("^rs", "rs123"));
        assert!(!is_match("^rs", "xrs123"));
        assert!(is_match("a+b", "a+b"));
        assert!(!is_match("a+b", "aab"));
        assert!(is_match("a\\+b", "aab"));
        assert!(is_match("1|0", "1|0"));
        assert!(!is_match("1|0", "1"));
        assert!(is_match("\\(ab\\)*c", "ababc"));
        assert!(is_match("*a", "*a"));
        assert!(is_match("a\\{2\\}", "aa"));
        assert!(!is_match("a\\{2\\}", "a"));
        assert!(is_match("[[:digit:]]", "x1"));
        assert!(is_match("[\\]", "\\"));
        assert!(is_match("a$b", "a$b"));
        assert!(is_match("INDEL", "SNP,INDEL"));
        assert!(is_match("<", "<DEL>"));
        assert!(is_match("[<>]", "<"));
    }

    #[test]
    fn back_references() {
        // As glibc's regexec matches them.
        let cases = [
            ("^\\(rs[0-9]*\\)\\1$", "rs12rs12", true),
            ("^\\(rs[0-9]*\\)\\1$", "rs12rs13", false),
            ("^\\(a*\\)\\1$", "aaaa", true),
            ("^\\(a*\\)\\1$", "aaa", false),
            ("^\\(.*\\)\\1$", "abab", true),
            ("x\\(a*\\)y\\1z", "xaayaaz", true),
            ("\\(a\\)\\1*b", "ab", true),
            ("\\(a\\)\\1\\+b", "ab", false),
            ("\\(a\\)\\1\\+b", "aaab", true),
            // One digit only: `\\10` is `\\1` and `0`.
            ("\\(a\\)\\10", "a10", false),
            ("\\(a\\)\\10", "aa0", true),
            // Nested groups, numbered as they open.
            ("^\\(\\(a\\)b\\)\\2\\1$", "abaab", true),
            ("^\\(\\(a\\)b\\)*c\\2$", "ababca", true),
            ("^\\(\\(a\\)b\\)*c\\2$", "c", false),
            ("^\\(\\(a\\)\\|b\\)*\\2$", "aba", true),
            ("^\\(\\(a\\)\\|b\\)*\\2$", "abb", false),
            // A group repeated no times keeps its number, as in glibc.
            ("\\(a\\)\\{0\\}\\(b\\)\\2", "bb", true),
            ("\\(a\\)\\{0\\}\\(b\\)\\1", "bb", false),
            ("\\(a\\)\\{0,0\\}\\(b\\)\\2", "bb", true),
            // A group that has not matched matches nothing.
            ("\\(a\\)*b\\1", "b", false),
            ("^\\(\\(a\\)\\|b\\)\\2$", "b", false),
            ("^\\(\\(a\\)\\|b\\)\\2$", "aa", true),
            // A group that matched the empty string, also in its last repetition.
            ("\\(\\)\\1", "x", true),
            ("^\\(\\)\\1$", "", true),
            ("^\\(a*\\)*b\\1$", "aab", true),
            ("^\\(a*\\)*b\\1$", "aaba", true),
            // `\\(` inside a bracket expression is two literal characters.
            ("^\\([\\(]\\)\\1$", "((", true),
            ("^\\([\\(]\\)\\1$", "\\\\", true),
            ("^\\([\\(]\\)\\1$", "(\\", false),
            ("\\<\\(a\\)\\1\\>", "x aa y", true),
            ("\\<\\(a\\)\\1\\>", "xaa y", false),
        ];
        for (p, s, want) in cases {
            assert_eq!(is_match(p, s), want, "{p} on {s:?}");
        }
        let icase = |p: &str, s: &str| {
            compile(p.as_bytes(), Case::Insensitive)
                .expect("the pattern compiles")
                .is_match(s.as_bytes())
        };
        assert!(icase("^\\(A\\)\\1$", "Aa"));
        assert!(icase("^\\(rs[0-9]\\)\\1$", "RS1rs1"));
        assert!(!is_match("^\\(A\\)\\1$", "Aa"));
    }

    #[test]
    fn back_references_on_long_texts() {
        // The matcher's stack is on the heap, so repetitions over 100 kb do not overflow the
        // thread's stack.
        let a = "a".repeat(100_000);
        let ab = format!("{a}b");
        assert!(is_match("\\(a\\)\\1.*b$", &ab));
        assert!(!is_match("\\(a\\)\\1a*$", &ab));
        assert!(is_match("^\\(a*\\)\\1b", &ab));
        assert!(!is_match("^\\(a*\\)\\1b", &format!("a{ab}")));
        assert!(is_match("\\(a\\)\\1\\(\\(\\(.\\)\\)\\)*b$", &ab));
        let sv = |last: char| format!("A{}{last}", "C".repeat(100_000));
        assert!(is_match("^\\(.\\).*\\1$", &sv('A')));
        assert!(!is_match("^\\(.\\).*\\1$", &sv('G')));
        assert!(is_match("\\(.\\).*\\1G$", &sv('G')));
        assert!(!is_match("\\(A\\).*\\1G$", &sv('G')));
    }

    #[test]
    fn back_references_in_nested_repetitions() {
        // Each state is tried once, so the iterations of the inner repetition are not tried
        // over every partition of the text (2^40 here).
        let a = "a".repeat(40);
        assert!(!is_match("^\\(a*\\)*\\1x$", &format!("{a}bx")));
        assert!(is_match("^\\(a*\\)*\\1x$", &format!("{a}x")));
        assert!(!is_match("^\\(\\(a\\)*\\|b\\)*\\2x$", &format!("{a}bbx")));
        assert!(!is_match("^\\(\\(a*\\)*\\)*\\1\\2x", &format!("{a}bx")));
        let a = "a".repeat(300);
        assert!(!is_match("^\\(a*\\)*\\1x$", &format!("{a}bx")));
        assert!(is_match("^\\(a*\\)*\\1b\\1x$", &format!("{a}bx")));
    }

    #[test]
    fn back_references_to_groups_closed_before() {
        // glibc refuses a reference to a group not closed before it in the same alternative.
        for p in [
            "\\1",
            "\\(a\\1\\)",
            "\\1\\(a\\)",
            "\\(a\\)\\|b\\1",
            "[\\(]\\1",
            "\\(a\\)\\2",
        ] {
            assert!(compile(p.as_bytes(), Case::Sensitive).is_err(), "{p}");
        }
        assert!(is_match("\\(\\(a\\)\\|b\\)\\2", "aa"));
        assert!(is_match("\\(a\\)\\|\\(b\\)\\2", "bb"));
    }

    #[test]
    fn only_back_references_backtrack() {
        let automaton = |p: &str| {
            let bre = compile(p.as_bytes(), Case::Sensitive).expect("the pattern compiles");
            matches!(bre.0, Matcher::Automaton(_))
        };
        assert!(automaton("^\\(rs[0-9]*\\)x$"));
        assert!(automaton("[\\1]"));
        assert!(!automaton("^\\(rs[0-9]*\\)\\1$"));
    }

    #[test]
    fn repetition_operators_after_anchors_are_literal() {
        // As in glibc, whose parser returns after an anchor: a `*`, `\+`, or `\?` after one
        // starts the next expression, as a literal character, and a `^` after one is literal.
        let cases = [
            ("b\\b*a", "b*a", true),
            ("b\\b*a", "ba", false),
            ("a\\>*", "a*", true),
            ("a\\>*", "a", false),
            ("\\`*a", "*a", true),
            ("\\`*a", "a", false),
            ("\\<*a", "a", false),
            ("a\\b\\+", "a+", true),
            ("a\\b\\+", "a", false),
            ("a\\'\\?", "a", false),
            ("^^a", "^a", true),
            ("^^a", "a", false),
            ("\\(^^a\\)", "^a", true),
            ("x\\|^^a", "a", false),
            ("\\b^*a", "a", true),
            ("\\b^*a", "^^a", true),
            ("\\b**a", "a", true),
            ("\\b*\\{2\\}a", "b**a", true),
            ("\\b*\\{2\\}a", "b*a", false),
            // Also with back-references.
            ("b\\b*\\(a\\)\\1", "b*aa", true),
            ("b\\b*\\(a\\)\\1", "baa", false),
            ("\\(a\\)\\1\\>\\+", "aa+", true),
            ("\\(a\\)\\1\\>\\+", "aa", false),
        ];
        for (p, s, want) in cases {
            assert_eq!(is_match(p, s), want, "{p} on {s:?}");
        }
        assert!(
            compile(b"b\\b*a", Case::Insensitive)
                .expect("the pattern compiles")
                .is_match(b"B*A")
        );
        // glibc refuses a counted repetition of an anchor.
        for p in [
            "^\\>\\{1\\}a",
            "\\b\\{1\\}a",
            "\\(a\\)\\1\\<\\{2\\}",
            "\\'\\{1,\\}",
        ] {
            assert!(compile(p.as_bytes(), Case::Sensitive).is_err(), "{p}");
        }
    }

    #[test]
    fn repetitions_of_repetitions() {
        // glibc applies a `\+` or `\?` after a repetition to the repeated expression.
        let cases = [
            ("a\\{2\\}\\?c", "c", true),
            ("a\\+\\?c", "c", true),
            ("a\\?\\+c", "aac", true),
            ("a\\+\\+c", "c", false),
            ("a\\{1,2\\}\\?\\+c", "c", true),
            ("\\(ab\\)\\{2\\}\\?c", "abc", true),
            ("[ab]\\+\\?c", "c", true),
            ("a*\\+b", "b", true),
            ("**a", "a", true),
            ("x\\(a\\)\\{2\\}\\?\\1", "xaaa", true),
            ("x\\(a\\)\\{2\\}\\?\\1", "xa", false),
            ("\\(a\\)\\(b\\)\\+\\?\\2c", "abbc", true),
            ("\\(a\\)\\(b\\)\\+\\?\\2c", "ac", false),
        ];
        for (p, s, want) in cases {
            assert_eq!(is_match(p, s), want, "{p} on {s:?}");
        }
        // It refuses a `*` or `\{` after one.
        for p in [
            "ab**",
            "a\\+*",
            "a\\?*",
            "a\\{1\\}*",
            "a*\\{2\\}",
            "a\\+\\{2\\}",
            "***a",
            "\\b***a",
            "\\(a\\)\\1*\\{2\\}",
        ] {
            assert!(compile(p.as_bytes(), Case::Sensitive).is_err(), "{p}");
        }
    }

    #[test]
    fn intervals() {
        // glibc reads a missing lower bound as 0.
        let cases = [
            ("^a\\{,2\\}c$", "aac", true),
            ("^a\\{,2\\}c$", "aaac", false),
            ("^a\\{,02\\}c$", "c", true),
            ("^a\\{,0\\}c$", "ac", false),
            ("^a\\{,\\}c$", "aaaac", true),
            ("^a\\{,\\}c$", "c", true),
            ("^a\\{00,2\\}c$", "c", true),
            ("^[ab]\\{,2\\}c$", "bac", true),
            ("^\\(a\\{,2\\}\\)\\1c$", "aaaac", true),
            ("^\\(a\\{,2\\}\\)\\1c$", "aaac", false),
            ("^a\\{32767\\}c$", "ac", false),
            // `\\,` and `\\0` are glibc's characters `,` and `0`.
            ("^a\\{1\\,2\\}c$", "aac", true),
            ("^a\\{\\02\\}c$", "aac", true),
            ("^a\\{,\\02\\}c$", "c", true),
            ("^\\(a\\)\\{1\\,\\}\\1c$", "aac", true),
        ];
        for (p, s, want) in cases {
            assert_eq!(is_match(p, s), want, "{p} on {s:?}");
        }
        // It refuses a count over RE_DUP_MAX, a lower bound over the upper one, and any other
        // content.
        for p in [
            "a\\{32768\\}",
            "a\\{0,32768\\}",
            "a\\{,32768\\}",
            "a\\{32768,\\}",
            "a\\{99999999999\\}",
            "a\\{2,1\\}",
            "a\\{\\}",
            "a\\{1,2,3\\}",
            "a\\{,1,\\}",
            "a\\{1a\\}",
            "a\\{ 1\\}",
            "a\\{1",
            "a\\{1,\\,\\}",
            "a\\{1\\2\\}",
        ] {
            assert!(compile(p.as_bytes(), Case::Sensitive).is_err(), "{p}");
        }
    }

    #[test]
    fn case_insensitive() {
        assert!(
            compile(b"pacbio", Case::Insensitive)
                .expect("the pattern compiles")
                .is_match(b"PacBio")
        );
    }
}
