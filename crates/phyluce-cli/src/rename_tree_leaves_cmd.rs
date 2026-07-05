//! CLI wiring for `phyluce genetrees rename-tree-leaves`, mirroring
//! `phyluce_genetrees_rename_tree_leaves`.
//!
//! `--reroot` is not yet implemented (it requires physically restructuring
//! the parsed tree, which none of the other genetree algorithms need).

use std::collections::HashMap;
use std::path::Path;

use phyluce_genetrees::newick::{parse_all, rename_leaves, write};

use crate::conf::parse_ini;

pub fn run(
    input: &Path,
    config: &Path,
    output: &Path,
    section: &str,
    order: &str,
    reroot: Option<&str>,
) -> anyhow::Result<()> {
    anyhow::ensure!(reroot.is_none(), "--reroot is not yet implemented");

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
        out.push_str(&write(&renamed));
        out.push('\n');
    }
    std::fs::write(output, out)?;
    Ok(())
}
