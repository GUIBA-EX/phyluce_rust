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

use phyluce_genetrees::newick::{bipartitions_polarized_by, parse, strip_branch_lengths, write};

pub fn run(
    trees_dir: &Path,
    locus_support_output: &Path,
    root: &str,
    extension: &str,
    exclude: &[String],
) -> anyhow::Result<()> {
    println!("creating tree objects");
    let mut entries: Vec<std::path::PathBuf> = std::fs::read_dir(trees_dir)?
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
        let text = std::fs::read_to_string(&tree_path)?;
        loci.push(LocusTree {
            gene_name: name,
            tree_text: text,
        });
    }

    println!("creating treelist");
    // key: sorted bipartition set (as a stable string); value: (count,
    // representative newick, loci names)
    let mut groups: HashMap<Vec<Vec<String>>, (usize, String, Vec<String>)> = HashMap::new();

    for locus in &loci {
        let tree = parse(&locus.tree_text)?;
        let stripped = strip_branch_lengths(&tree);
        let key: Vec<Vec<String>> = bipartitions_polarized_by(&stripped, root)
            .into_iter()
            .collect();
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
        println!("{count}\t{newick}");
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
    std::fs::write(locus_support_output, out)?;
    Ok(())
}
