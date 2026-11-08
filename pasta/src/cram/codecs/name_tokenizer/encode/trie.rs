//! A trie of read names that finds, for each new name, the last earlier name through a node of its
//! path: htscodecs' `build_trie` and `search_trie`, with each run of nodes that no name branches
//! from or ends in merged into one.

const NONE: u32 = u32::MAX;
const ROOT: usize = 0;

pub(super) struct Trie<'a> {
    nodes: Vec<Node<'a>>,
}

struct Node<'a> {
    /// The bytes on the edge into this node. Every name through the node passes through all of
    /// them, so the nodes they stand for have the same last visitor as this one.
    label: &'a [u8],
    first_child: u32,
    next_sibling: u32,
    /// The last name through this node.
    last: u32,
}

/// The last names through two nodes of a name's path, before the name itself.
#[derive(Debug, Default, PartialEq)]
pub(super) struct Visit {
    /// The last name through the node of the whole name: an equal name or one that the name
    /// begins, else the name itself. `None` for an empty name.
    pub(super) full: Option<u32>,
    /// The last name through the node `prefix_len` bytes in, as for [`Self::full`]. `None` if the
    /// name is shorter.
    pub(super) prefix: Option<u32>,
}

impl Visit {
    /// Records that the nodes at depths `start + 1..=end` of `name` were last visited by `last`.
    fn reach(&mut self, start: usize, end: usize, name: &[u8], prefix_len: usize, last: u32) {
        if start < prefix_len && prefix_len <= end {
            self.prefix = Some(last);
        }

        if end == name.len() {
            self.full = Some(last);
        }
    }
}

impl<'a> Trie<'a> {
    pub(super) fn new() -> Self {
        Self {
            nodes: vec![Node::new(&[], NONE)],
        }
    }

    /// Adds name `n`, returning the last names before it through its nodes.
    pub(super) fn visit(&mut self, name: &'a [u8], n: u32, prefix_len: usize) -> Visit {
        let mut visit = Visit::default();
        let mut node = ROOT;
        let mut depth = 0;

        while let Some(&b) = name.get(depth) {
            let mut prev_sibling = None;
            let mut child = self.nodes[node].first_child;

            while child != NONE && self.nodes[child as usize].label[0] != b {
                prev_sibling = Some(child as usize);
                child = self.nodes[child as usize].next_sibling;
            }

            let rest = &name[depth..];

            if child == NONE {
                let leaf = self.push(rest, n);
                self.link(node, prev_sibling, leaf);
                visit.reach(depth, name.len(), name, prefix_len, n);
                return visit;
            }

            let child = child as usize;
            let label = self.nodes[child].label;
            let len = common_prefix_len(label, rest);
            visit.reach(depth, depth + len, name, prefix_len, self.nodes[child].last);

            if len == label.len() {
                self.nodes[child].last = n;
                node = child;
                depth += len;
                continue;
            }

            // The name ends or branches off inside the edge, which is split there.
            let mid = self.push(&label[..len], n);
            self.nodes[mid].first_child = child as u32;
            self.nodes[mid].next_sibling = self.nodes[child].next_sibling;
            self.link(node, prev_sibling, mid);

            self.nodes[child].label = &label[len..];
            self.nodes[child].next_sibling = NONE;

            depth += len;

            if depth < name.len() {
                let leaf = self.push(&name[depth..], n);
                self.nodes[child].next_sibling = leaf as u32;
                visit.reach(depth, name.len(), name, prefix_len, n);
            }

            return visit;
        }

        visit
    }

    fn push(&mut self, label: &'a [u8], n: u32) -> usize {
        self.nodes.push(Node::new(label, n));
        self.nodes.len() - 1
    }

    /// Puts `node` in place of `parent`'s first child, or of the sibling after `prev_sibling`.
    fn link(&mut self, parent: usize, prev_sibling: Option<usize>, node: usize) {
        let node = node as u32;

        match prev_sibling {
            Some(i) => self.nodes[i].next_sibling = node,
            None => self.nodes[parent].first_child = node,
        }
    }
}

impl<'a> Node<'a> {
    fn new(label: &'a [u8], last: u32) -> Self {
        Self {
            label,
            first_child: NONE,
            next_sibling: NONE,
            last,
        }
    }
}

fn common_prefix_len(a: &[u8], b: &[u8]) -> usize {
    a.iter().zip(b).take_while(|(x, y)| x == y).count()
}

#[cfg(test)]
mod tests {
    use super::*;

    /// The last name before `n` that starts with `name[..depth]`, or `n` if none does.
    fn last_through(names: &[&[u8]], n: usize, depth: usize) -> Option<u32> {
        let name = names[n];

        (depth > 0 && depth <= name.len()).then(|| {
            let prefix = &name[..depth];
            let last = names[..n]
                .iter()
                .rposition(|other| other.starts_with(prefix));
            last.unwrap_or(n) as u32
        })
    }

    #[test]
    fn test_visit() {
        let mut names: Vec<Vec<u8>> = [
            "A:1:2:3:4",
            "A:1:2:3:4",
            "A:1:2:3:45",
            "A:1:2:3:4",
            "A:1:2:3",
            "B",
            "",
            "A:1:2:3:45",
            "A:1:9",
            "AB",
        ]
        .iter()
        .map(|s| s.as_bytes().to_vec())
        .collect();

        // Names that share prefixes of every length.
        let mut x = 0x9e37_79b9_u32;

        for _ in 0..2000 {
            x ^= x << 13;
            x ^= x >> 17;
            x ^= x << 5;

            let len = (x % 9) as usize;
            names.push(
                (0..len)
                    .map(|i| b"ab:"[((x >> (i * 2)) % 3) as usize])
                    .collect(),
            );
        }

        let names: Vec<&[u8]> = names.iter().map(|name| &name[..]).collect();

        for prefix_len in [1, 2, 4, 7, usize::MAX] {
            let mut trie = Trie::new();

            for (n, &name) in names.iter().enumerate() {
                let expected = Visit {
                    full: last_through(&names, n, name.len()),
                    prefix: last_through(&names, n, prefix_len),
                };

                let actual = trie.visit(name, n as u32, prefix_len);
                assert_eq!(
                    actual, expected,
                    "name {n} {name:?}, prefix length {prefix_len}"
                );
            }
        }
    }
}
