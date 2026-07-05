//! A minimal Newick tree parser/writer, sufficient for the genetrees
//! commands: leaf renaming, bootstrap-support-label extraction, and
//! rooted-topology comparison. Mirrors DendroPy's `preserve_underscores=True`
//! (underscores in unquoted labels are kept literal, not converted to
//! spaces per the strict Newick spec).

use std::collections::{BTreeSet, HashMap};

#[derive(Debug, Clone, PartialEq)]
pub struct Node {
    /// Leaf name, or an internal node's label (commonly a bootstrap
    /// support value written by RAxML).
    pub label: Option<String>,
    pub branch_length: Option<f64>,
    pub children: Vec<Node>,
}

impl Node {
    pub fn is_leaf(&self) -> bool {
        self.children.is_empty()
    }
}

#[derive(Debug, thiserror::Error)]
pub enum NewickError {
    #[error("unexpected end of input while parsing Newick tree")]
    UnexpectedEof,
    #[error("expected '{expected}' at byte {pos}, found {found:?}")]
    Expected {
        expected: char,
        pos: usize,
        found: Option<char>,
    },
}

struct Parser<'a> {
    chars: Vec<char>,
    pos: usize,
    _marker: std::marker::PhantomData<&'a ()>,
}

impl<'a> Parser<'a> {
    fn new(text: &'a str) -> Self {
        Parser {
            chars: text.chars().collect(),
            pos: 0,
            _marker: std::marker::PhantomData,
        }
    }

    fn peek(&self) -> Option<char> {
        self.chars.get(self.pos).copied()
    }

    fn advance(&mut self) -> Option<char> {
        let c = self.peek();
        if c.is_some() {
            self.pos += 1;
        }
        c
    }

    fn skip_ws(&mut self) {
        while matches!(self.peek(), Some(c) if c.is_whitespace()) {
            self.pos += 1;
        }
    }

    fn expect(&mut self, expected: char) -> Result<(), NewickError> {
        self.skip_ws();
        match self.advance() {
            Some(c) if c == expected => Ok(()),
            found => Err(NewickError::Expected {
                expected,
                pos: self.pos,
                found,
            }),
        }
    }

    fn parse_label(&mut self) -> Option<String> {
        self.skip_ws();
        if self.peek() == Some('\'') {
            self.advance();
            let mut s = String::new();
            while let Some(c) = self.advance() {
                if c == '\'' {
                    if self.peek() == Some('\'') {
                        s.push('\'');
                        self.advance();
                        continue;
                    }
                    break;
                }
                s.push(c);
            }
            return Some(s);
        }
        let mut s = String::new();
        while let Some(c) = self.peek() {
            if "(),:;".contains(c) || c.is_whitespace() {
                break;
            }
            s.push(c);
            self.pos += 1;
        }
        if s.is_empty() {
            None
        } else {
            Some(s)
        }
    }

    fn parse_branch_length(&mut self) -> Option<f64> {
        self.skip_ws();
        if self.peek() != Some(':') {
            return None;
        }
        self.advance();
        self.skip_ws();
        let mut s = String::new();
        while let Some(c) = self.peek() {
            if c.is_ascii_digit() || c == '.' || c == '-' || c == '+' || c == 'e' || c == 'E' {
                s.push(c);
                self.pos += 1;
            } else {
                break;
            }
        }
        s.parse::<f64>().ok()
    }

    fn parse_subtree(&mut self) -> Result<Node, NewickError> {
        self.skip_ws();
        let children = if self.peek() == Some('(') {
            self.advance();
            let mut kids = vec![self.parse_subtree()?];
            self.skip_ws();
            while self.peek() == Some(',') {
                self.advance();
                kids.push(self.parse_subtree()?);
                self.skip_ws();
            }
            self.expect(')')?;
            kids
        } else {
            Vec::new()
        };
        let label = self.parse_label();
        let branch_length = self.parse_branch_length();
        Ok(Node {
            label,
            branch_length,
            children,
        })
    }
}

/// Parse a single Newick tree (terminated by `;`, which is consumed if
/// present).
pub fn parse(text: &str) -> Result<Node, NewickError> {
    let mut p = Parser::new(text.trim());
    let node = p.parse_subtree()?;
    p.skip_ws();
    if p.peek() == Some(';') {
        p.advance();
    }
    Ok(node)
}

/// Parse every semicolon-terminated tree in a multi-line/multi-tree file.
pub fn parse_all(text: &str) -> Result<Vec<Node>, NewickError> {
    text.split(';')
        .map(str::trim)
        .filter(|s| !s.is_empty())
        .map(|s| parse(&format!("{s};")))
        .collect()
}

fn write_node(node: &Node, out: &mut String) {
    if !node.children.is_empty() {
        out.push('(');
        for (i, child) in node.children.iter().enumerate() {
            if i > 0 {
                out.push(',');
            }
            write_node(child, out);
        }
        out.push(')');
    }
    if let Some(label) = &node.label {
        out.push_str(label);
    }
    if let Some(bl) = node.branch_length {
        out.push(':');
        out.push_str(&format!("{bl}"));
    }
}

pub fn write(node: &Node) -> String {
    let mut s = String::new();
    write_node(node, &mut s);
    s.push(';');
    s
}

pub fn strip_branch_lengths(node: &Node) -> Node {
    Node {
        label: node.label.clone(),
        branch_length: None,
        children: node.children.iter().map(strip_branch_lengths).collect(),
    }
}

pub fn leaves(node: &Node) -> BTreeSet<String> {
    let mut out = BTreeSet::new();
    collect_leaves(node, &mut out);
    out
}

fn collect_leaves(node: &Node, out: &mut BTreeSet<String>) {
    if node.is_leaf() {
        if let Some(l) = &node.label {
            out.insert(l.clone());
        }
    } else {
        for c in &node.children {
            collect_leaves(c, out);
        }
    }
}

pub fn rename_leaves(node: &Node, mapping: &HashMap<String, String>) -> Node {
    if node.is_leaf() {
        let new_label = node.label.as_ref().and_then(|l| {
            mapping
                .get(l)
                .or_else(|| mapping.get(&l.replace(' ', "_")))
                .cloned()
        });
        Node {
            label: new_label.or_else(|| node.label.clone()),
            branch_length: node.branch_length,
            children: Vec::new(),
        }
    } else {
        Node {
            label: node.label.clone(),
            branch_length: node.branch_length,
            children: node
                .children
                .iter()
                .map(|c| rename_leaves(c, mapping))
                .collect(),
        }
    }
}

/// Every internal node's descendant leaf set, side-normalized so it never
/// contains `outgroup` -- a rooting-invariant representation sufficient to
/// compare two trees' topologies once both are considered "rooted at (or
/// polarized by) the same outgroup taxon". Mirrors comparing
/// `tree.reroot_at_edge(...)`-then-`symmetric_difference` in spirit,
/// without needing to physically mutate the tree.
pub fn bipartitions_polarized_by(node: &Node, outgroup: &str) -> BTreeSet<Vec<String>> {
    let all_leaves = leaves(node);
    let mut out = BTreeSet::new();
    collect_bipartitions(node, &all_leaves, outgroup, &mut out);
    out
}

fn collect_bipartitions(
    node: &Node,
    all_leaves: &BTreeSet<String>,
    outgroup: &str,
    out: &mut BTreeSet<Vec<String>>,
) {
    if node.is_leaf() {
        return;
    }
    let side = leaves(node);
    // non-trivial bipartitions only (every leaf, or a single leaf, carries
    // no topological information).
    if side.len() > 1 && side.len() < all_leaves.len() {
        let normalized: Vec<String> = if side.contains(outgroup) {
            all_leaves.difference(&side).cloned().collect()
        } else {
            side.iter().cloned().collect()
        };
        out.insert(normalized);
    }
    for c in &node.children {
        collect_bipartitions(c, all_leaves, outgroup, out);
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn parses_and_writes_simple_tree() {
        let tree = parse("(a,b,(c,d));").unwrap();
        assert_eq!(
            leaves(&tree),
            BTreeSet::from(["a", "b", "c", "d"].map(String::from))
        );
        assert_eq!(write(&tree), "(a,b,(c,d));");
    }

    #[test]
    fn parses_labels_and_branch_lengths() {
        let tree = parse("(a:0.1,b:0.2)90:0.05;").unwrap();
        assert_eq!(tree.label.as_deref(), Some("90"));
        assert_eq!(tree.branch_length, Some(0.05));
        assert_eq!(tree.children[0].branch_length, Some(0.1));
    }

    #[test]
    fn renames_leaves() {
        let tree = parse("(a,b,(c,d));").unwrap();
        let map: HashMap<String, String> = [("a", "alpha"), ("b", "beta")]
            .into_iter()
            .map(|(k, v)| (k.to_string(), v.to_string()))
            .collect();
        let renamed = rename_leaves(&tree, &map);
        assert_eq!(
            leaves(&renamed),
            BTreeSet::from(["alpha", "beta", "c", "d"].map(String::from))
        );
    }

    #[test]
    fn same_topology_regardless_of_input_rooting() {
        let t1 = parse("((a,b),c,d);").unwrap();
        let t2 = parse("((c,d),a,b);").unwrap();
        assert_eq!(
            bipartitions_polarized_by(&t1, "d"),
            bipartitions_polarized_by(&t2, "d")
        );
    }

    #[test]
    fn different_topology_detected() {
        let t1 = parse("((a,b),c,d);").unwrap();
        let t2 = parse("((a,c),b,d);").unwrap();
        assert_ne!(
            bipartitions_polarized_by(&t1, "d"),
            bipartitions_polarized_by(&t2, "d")
        );
    }
}
