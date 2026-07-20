//! CLI wiring for `phyluce genetrees rename-tree-leaves`, mirroring
//! `phyluce_genetrees_rename_tree_leaves`.

use std::collections::HashMap;
use std::path::Path;

use phyluce_genetrees::newick::{parse_all, rename_leaves, reroot_at_leaf_parent, write};

use crate::conf::parse_ini;

pub fn run(
    input: &Path,
    config: &Path,
    output: &Path,
    section: &str,
    order: &str,
    reroot: Option<&str>,
) -> anyhow::Result<()> {
    let conf_text = std::fs::read_to_string(config)?;
    let sections = parse_ini(&conf_text);
    let entries = sections
        .get(section)
        .ok_or_else(|| anyhow::anyhow!("no [{section}] section in --config"))?;

    let names: HashMap<String, String> = match order {
        "right:left" => entries
            .iter()
            .map(|(k, v)| (v.replace('-', "_"), k.clone()))
            .collect(),
        _ => entries
            .iter()
            .map(|(k, v)| (k.replace('-', "_"), v.clone()))
            .collect(),
    };

    let tree_text = std::fs::read_to_string(input)?;
    let trees = parse_all(&tree_text)?;
    let mut out = String::new();
    for tree in &trees {
        let renamed = rename_leaves(tree, &names);
        let final_tree = match reroot {
            Some(taxon) => reroot_at_leaf_parent(&renamed, taxon).ok_or_else(|| {
                anyhow::anyhow!("--reroot taxon {taxon:?} not found among tree leaves")
            })?,
            None => renamed,
        };
        out.push_str(&write(&final_tree));
        out.push('\n');
    }
    std::fs::write(output, out)?;
    Ok(())
}
