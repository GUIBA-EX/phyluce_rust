//! CLI wiring for `phyluce genetrees get-tree-counts`, mirroring
//! `phyluce_genetrees_get_tree_counts`.
//!
//! Topology grouping uses a rooting-invariant bipartition comparison
//! (`phyluce_genetrees::newick::bipartitions_polarized_by`) rather than
//! physically rerooting each tree at `--root` the way the Python original
//! does, so the printed representative Newick strings keep their
//! as-parsed rooting/branch order instead of being re-rooted at the
//! chosen outgroup. Grouping/counting results are equivalent; only the
//! cosmetic rooting of the printed tree strings differs.

use std::collections::HashMap;
use std::path::Path;

use anyhow::Context;
use phyluce_genetrees::newick::{
    bipartitions_polarized_by, leaves, parse, strip_branch_lengths, write, Node,
};

type TopologyKey = (std::collections::BTreeSet<String>, Vec<Vec<String>>);
type TopologyGroup = (usize, String, Vec<String>);

fn topology_key(tree: &Node, root: &str, locus: &str) -> anyhow::Result<TopologyKey> {
    let stripped = strip_branch_lengths(tree);
    let leaf_set = leaves(&stripped);
    anyhow::ensure!(
        leaf_set.contains(root),
        "outgroup {root:?} is absent from tree for locus {locus}"
    );
    Ok((
        leaf_set,
        bipartitions_polarized_by(&stripped, root)
            .into_iter()
            .collect(),
    ))
}

pub fn run(
    trees_dir: &Path,
    locus_support_output: &Path,
    root: &str,
    extension: &str,
    exclude: &[String],
) -> anyhow::Result<()> {
    crate::cli_info!("creating tree objects");
    let mut entries: Vec<std::path::PathBuf> = std::fs::read_dir(trees_dir)
        .with_context(|| format!("reading trees directory {}", trees_dir.display()))?
        .filter_map(|e| e.ok())
        .map(|e| e.path())
        .filter(|p| p.is_dir())
        .collect();
    entries.sort();

    struct LocusTree {
        gene_name: String,
        tree_text: String,
    }
    let mut loci = Vec::new();
    for dir in &entries {
        let name = dir
            .file_name()
            .and_then(|n| n.to_str())
            .unwrap_or("")
            .to_string();
        if exclude.contains(&name) {
            continue;
        }
        let tree_path = dir.join(format!("RAxML_bestTree.{extension}"));
        if !tree_path.is_file() {
            continue;
        }
        let text = std::fs::read_to_string(&tree_path)
            .with_context(|| format!("reading tree file {}", tree_path.display()))?;
        loci.push(LocusTree {
            gene_name: name,
            tree_text: text,
        });
    }

    crate::cli_info!("creating treelist");
    // The leaf set is part of the key so incomplete trees cannot be grouped
    // solely because they happen to have the same remaining bipartitions.
    let mut groups: HashMap<TopologyKey, TopologyGroup> = HashMap::new();

    for locus in &loci {
        let tree = parse(&locus.tree_text)?;
        let stripped = strip_branch_lengths(&tree);
        let key = topology_key(&tree, root, &locus.gene_name)?;
        let entry = groups
            .entry(key)
            .or_insert_with(|| (0, write(&stripped), Vec::new()));
        entry.0 += 1;
        entry.2.push(locus.gene_name.clone());
    }

    let mut by_count: Vec<(usize, String)> =
        groups.values().map(|(c, s, _)| (*c, s.clone())).collect();
    by_count.sort_by_key(|b| std::cmp::Reverse(b.0));
    for (count, newick) in &by_count {
        crate::cli_info!("{count}\t{newick}");
    }

    let mut by_loci: Vec<(String, Vec<String>)> = groups
        .into_values()
        .map(|(_, newick, loci)| (newick, loci))
        .collect();
    by_loci.sort_by_key(|b| std::cmp::Reverse(b.1.len()));

    let mut out = String::new();
    for (i, (newick, loci)) in by_loci.iter().enumerate() {
        out.push_str(&format!(
            "# {}\n# {}\n[{}th most numerous]\n{}\n\n",
            newick,
            loci.len(),
            i + 1,
            loci.join("\n")
        ));
    }
    std::fs::write(locus_support_output, out).with_context(|| {
        format!(
            "writing locus support output to {}",
            locus_support_output.display()
        )
    })?;
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn topology_key_requires_the_requested_outgroup() {
        let tree = parse("(a,b,c);").unwrap();
        assert!(topology_key(&tree, "missing", "locus-1").is_err());
    }

    #[test]
    fn topology_key_keeps_leaf_sets_distinct() {
        let first = parse("(a,b,c);").unwrap();
        let second = parse("(d,e,f);").unwrap();
        assert_ne!(
            topology_key(&first, "a", "locus-1").unwrap(),
            topology_key(&second, "d", "locus-2").unwrap()
        );
    }
}
