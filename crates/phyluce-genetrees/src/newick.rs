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
    #[error("unterminated quoted label at character {pos}")]
    UnterminatedQuotedLabel { pos: usize },
    #[error("invalid branch length {value:?} at character {pos}")]
    InvalidBranchLength { value: String, pos: usize },
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

    fn parse_label(&mut self) -> Result<Option<String>, NewickError> {
        self.skip_ws();
        if self.peek() == Some('\'') {
            let start = self.pos;
            self.advance();
            let mut s = String::new();
            while let Some(c) = self.advance() {
                if c == '\'' {
                    if self.peek() == Some('\'') {
                        s.push('\'');
                        self.advance();
                        continue;
                    }
                    return Ok(Some(s));
                }
                s.push(c);
            }
            return Err(NewickError::UnterminatedQuotedLabel { pos: start });
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
            Ok(None)
        } else {
            Ok(Some(s))
        }
    }

    fn parse_branch_length(&mut self) -> Result<Option<f64>, NewickError> {
        self.skip_ws();
        if self.peek() != Some(':') {
            return Ok(None);
        }
        self.advance();
        self.skip_ws();
        let start = self.pos;
        let mut s = String::new();
        while let Some(c) = self.peek() {
            if c.is_ascii_digit() || c == '.' || c == '-' || c == '+' || c == 'e' || c == 'E' {
                s.push(c);
                self.pos += 1;
            } else {
                break;
            }
        }
        let value = s
            .parse::<f64>()
            .map_err(|_| NewickError::InvalidBranchLength {
                value: s,
                pos: start,
            })?;
        if !value.is_finite() {
            return Err(NewickError::InvalidBranchLength {
                value: value.to_string(),
                pos: start,
            });
        }
        Ok(Some(value))
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
        let label = self.parse_label()?;
        let branch_length = self.parse_branch_length()?;
        Ok(Node {
            label,
            branch_length,
            children,
        })
    }
}

/// Parse exactly one semicolon-terminated Newick tree.
pub fn parse(text: &str) -> Result<Node, NewickError> {
    let mut p = Parser::new(text);
    let node = p.parse_subtree()?;
    p.expect(';')?;
    p.skip_ws();
    if p.peek().is_some() {
        return Err(NewickError::Expected {
            expected: ';',
            pos: p.pos,
            found: p.peek(),
        });
    }
    Ok(node)
}

/// Parse every semicolon-terminated tree in a multi-line/multi-tree file.
pub fn parse_all(text: &str) -> Result<Vec<Node>, NewickError> {
    let mut parser = Parser::new(text);
    let mut trees = Vec::new();
    loop {
        parser.skip_ws();
        if parser.peek().is_none() {
            return Ok(trees);
        }
        let tree = parser.parse_subtree()?;
        parser.expect(';')?;
        trees.push(tree);
    }
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
        write_label(label, out);
    }
    if let Some(bl) = node.branch_length {
        out.push(':');
        out.push_str(&format!("{bl}"));
    }
}

fn write_label(label: &str, out: &mut String) {
    let needs_quotes = label.is_empty()
        || label
            .chars()
            .any(|c| c.is_whitespace() || "(),:;[]'".contains(c));
    if needs_quotes {
        out.push('\'');
        out.push_str(&label.replace('\'', "''"));
        out.push('\'');
    } else {
        out.push_str(label);
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

/// Reroot the tree at the parent of the leaf named `leaf_label`, mirroring
/// DendroPy's `tree.reroot_at_node(leaf.parent_node)`: the leaf's parent
/// becomes the new root, keeping its own children, with the rest of the
/// tree grafted on as one extra child (the path back to the old root is
/// inverted edge-by-edge). Returns `None` if no leaf has that label.
pub fn reroot_at_leaf_parent(root: &Node, leaf_label: &str) -> Option<Node> {
    let path = find_leaf_parent_path(root, leaf_label)?;
    Some(reroot_at_path(root, &path))
}

/// Path of child indices from `node` down to the parent of the leaf named
/// `label` (empty if `node` itself is that parent).
fn find_leaf_parent_path(node: &Node, label: &str) -> Option<Vec<usize>> {
    if node
        .children
        .iter()
        .any(|c| c.is_leaf() && c.label.as_deref() == Some(label))
    {
        return Some(Vec::new());
    }
    for (i, child) in node.children.iter().enumerate() {
        if let Some(mut sub) = find_leaf_parent_path(child, label) {
            let mut full = vec![i];
            full.append(&mut sub);
            return Some(full);
        }
    }
    None
}

fn reroot_at_path(root: &Node, path: &[usize]) -> Node {
    if path.is_empty() {
        return root.clone();
    }
    // Walk root -> target, peeling each ancestor's branch toward the next
    // node into a "stripped" node (itself minus that child, carrying that
    // edge's length as its own -- it will hang off the inverted chain).
    let mut cur = root;
    let mut stripped: Vec<Node> = Vec::new();
    for &idx in path {
        let child = &cur.children[idx];
        let mut remaining = cur.children.clone();
        remaining.remove(idx);
        stripped.push(Node {
            label: cur.label.clone(),
            branch_length: child.branch_length,
            children: remaining,
        });
        cur = child;
    }
    // Fold the stripped ancestors from the root side inward, each becoming
    // a child of the next one down, producing the single subtree that
    // used to sit "above" the target.
    let mut iter = stripped.into_iter();
    let mut inverted = iter.next().expect("path is non-empty");
    for mut next in iter {
        next.children.push(inverted);
        inverted = next;
    }
    let mut new_root = cur.clone();
    new_root.branch_length = None;
    new_root.children.push(inverted);
    new_root
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
    fn reroot_is_noop_when_leaf_parent_is_already_root() {
        // root's direct children are [x, y, z]; y's parent is already the
        // root, so rerooting on y changes nothing.
        let tree = parse("(x,y,z);").unwrap();
        let rerooted = reroot_at_leaf_parent(&tree, "y").unwrap();
        assert_eq!(write(&rerooted), "(x,y,z);");
    }

    #[test]
    fn reroots_at_leaf_parent_nested() {
        // ((a,b),c,(d,e)); reroot on "d" -> parent of d is (d,e), which
        // becomes the new root: it keeps its own children (d, e) and gains
        // one extra child holding the inverted rest of the tree ((a,b)
        // and c, formerly siblings of (d,e) under the old root).
        let tree = parse("((a,b),c,(d,e));").unwrap();
        let rerooted = reroot_at_leaf_parent(&tree, "d").unwrap();
        assert_eq!(leaves(&rerooted), leaves(&tree));
        assert_eq!(rerooted.children.len(), 3);
        let labels: Vec<Option<&str>> = rerooted
            .children
            .iter()
            .map(|c| c.label.as_deref())
            .collect();
        assert!(labels.contains(&Some("d")));
        assert!(labels.contains(&Some("e")));
        let remainder = rerooted
            .children
            .iter()
            .find(|c| c.label.is_none())
            .expect("remainder child");
        assert_eq!(
            leaves(remainder),
            BTreeSet::from(["a", "b", "c"].map(String::from))
        );
    }

    #[test]
    fn reroot_missing_leaf_returns_none() {
        let tree = parse("(a,b,c);").unwrap();
        assert!(reroot_at_leaf_parent(&tree, "nope").is_none());
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

    #[test]
    fn rejects_trailing_input_and_unterminated_quotes() {
        assert!(parse("(a,b); trailing").is_err());
        assert!(parse("('a,b);").is_err());
        assert!(parse("(a,b):;").is_err());
        assert!(parse("(a,b):not-a-number;").is_err());
    }

    #[test]
    fn parses_multiple_trees_with_semicolons_in_quoted_labels() {
        let trees = parse_all("('a;b',c);(d,e);").unwrap();
        assert_eq!(trees.len(), 2);
        assert_eq!(
            leaves(&trees[0]),
            BTreeSet::from(["a;b".to_string(), "c".to_string()])
        );
    }

    #[test]
    fn writer_quotes_labels_that_require_escaping() {
        let tree = parse("('a b','c;d','e''f');").unwrap();
        let text = write(&tree);
        assert_eq!(text, "('a b','c;d','e''f');");
        assert_eq!(parse(&text).unwrap(), tree);
    }
}
