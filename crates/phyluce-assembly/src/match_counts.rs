//! Contig-count / matrix-membership pipeline mirroring
//! `phyluce_assembly_get_match_counts`'s non-`--optimize` path (complete and
//! incomplete matrix generation from `probe.matches.sqlite`).
//!
//! `--optimize`/`--random` (combinatorial best-subgroup search) isn't
//! implemented yet -- it's a rarely used, heavy (multiprocessing) feature
//! with no golden fixture; see docs/rust-rewrite-plan.md's phased command
//! priority list.

use std::collections::{HashMap, HashSet};
use std::path::Path;

use phyluce_io::sql::{ident, qualified_ident};
use rusqlite::Connection;

#[derive(Debug, thiserror::Error)]
pub enum MatchCountError {
    #[error("{0}")]
    Io(#[from] std::io::Error),
    #[error("{0}")]
    Sqlite(#[from] rusqlite::Error),
    #[error("no [{0}] section in taxon-list-config")]
    NoSuchGroup(String),
}

/// A `[section]`-per-line conf file, `allow_no_value=True`-style: each
/// non-empty, non-comment line under a section is an item name (an
/// optional `:value`/`=value` suffix is ignored). Mirrors
/// `configparser.RawConfigParser(allow_no_value=True)` as used for
/// `--taxon-list-config`.
pub fn read_taxon_list_config(path: &Path) -> std::io::Result<HashMap<String, Vec<String>>> {
    let text = std::fs::read_to_string(path)?;
    let mut sections: HashMap<String, Vec<String>> = HashMap::new();
    let mut current: Option<String> = None;
    for raw_line in text.lines() {
        let line = raw_line.trim();
        if line.is_empty() || line.starts_with('#') {
            continue;
        }
        if line.starts_with('[') && line.ends_with(']') {
            let name = line[1..line.len() - 1].trim().to_string();
            sections.entry(name.clone()).or_default();
            current = Some(name);
            continue;
        }
        if let Some(section) = &current {
            let key = line
                .split_once(':')
                .or_else(|| line.split_once('='))
                .map(|(k, _)| k.trim())
                .unwrap_or(line);
            sections
                .entry(section.clone())
                .or_default()
                .push(key.to_string());
        }
    }
    Ok(sections)
}

/// Mirrors `get_names_from_config`: taxon names for a group, `-` replaced
/// with `_`.
pub fn names_from_config(
    config: &HashMap<String, Vec<String>>,
    group: &str,
) -> Result<Vec<String>, MatchCountError> {
    match config.get(group) {
        Some(names) => Ok(names.iter().map(|n| n.replace('-', "_")).collect()),
        None if group == "Excludes" => Ok(Vec::new()),
        None => Err(MatchCountError::NoSuchGroup(group.to_string())),
    }
}

/// Mirrors `get_taxa_from_config`: names in `group`, minus any names listed
/// under `[Excludes]`.
pub fn taxa_from_config(
    config: &HashMap<String, Vec<String>>,
    group: &str,
) -> Result<Vec<String>, MatchCountError> {
    let organisms = names_from_config(config, group)?;
    let excludes: HashSet<String> = names_from_config(config, "Excludes")?.into_iter().collect();
    Ok(organisms
        .into_iter()
        .filter(|o| !excludes.contains(o))
        .collect())
}

/// Mirrors `get_uce_names`: every UCE locus name in `matches`.
pub fn uce_names(conn: &Connection) -> rusqlite::Result<HashSet<String>> {
    let mut stmt = conn.prepare("SELECT uce FROM matches")?;
    let rows = stmt.query_map([], |r| r.get::<_, String>(0))?;
    rows.collect()
}

/// Mirrors `remove_duplicates_from`: given the set of UCE loci an organism
/// matched, drop any locus whose `match_map` contig also matched another
/// locus for the same organism (a duplicate query/contig match).
fn remove_duplicates_from(
    conn: &Connection,
    organism: &str,
    matches: &HashSet<String>,
) -> rusqlite::Result<Vec<String>> {
    if matches.is_empty() {
        return Ok(Vec::new());
    }
    let (table, column) = if let Some(stripped) = organism.strip_suffix('*') {
        ("extended.match_map", stripped.to_string())
    } else {
        ("match_map", organism.to_string())
    };
    let in_list = vec!["?"; matches.len()].join(",");
    let column = ident(&column);
    let table = qualified_ident(table);
    let query = format!("SELECT uce, {column} FROM {table} WHERE uce in ({in_list})");
    let mut stmt = conn.prepare(&query)?;
    let rows: Vec<(String, String)> = stmt
        .query_map(rusqlite::params_from_iter(matches.iter()), |r| {
            Ok((r.get::<_, String>(0)?, r.get::<_, String>(1)?))
        })?
        .collect::<Result<_, _>>()?;

    let mut by_node: HashMap<String, Vec<String>> = HashMap::new();
    for (uce, contig_field) in rows {
        let node = contig_field.split('(').next().unwrap_or("").to_string();
        by_node.entry(node).or_default().push(uce);
    }
    Ok(by_node
        .into_values()
        .filter(|v| v.len() <= 1)
        .map(|v| v[0].clone())
        .collect())
}

/// Mirrors `get_all_matches_by_organism`.
pub fn matches_by_organism(
    conn: &Connection,
    organisms: &[String],
) -> rusqlite::Result<HashMap<String, HashSet<String>>> {
    let mut out = HashMap::new();
    for organism in organisms {
        let (table, column) = if let Some(stripped) = organism.strip_suffix('*') {
            ("extended.matches", stripped.to_string())
        } else {
            ("matches", organism.clone())
        };
        let table = qualified_ident(table);
        let column = ident(&column);
        let query = format!("SELECT uce FROM {table} WHERE {column} = 1");
        let mut stmt = conn.prepare(&query)?;
        let raw: HashSet<String> = stmt
            .query_map([], |r| r.get::<_, String>(0))?
            .collect::<Result<_, _>>()?;
        let deduped = remove_duplicates_from(conn, organism, &raw)?;
        out.insert(organism.clone(), deduped.into_iter().collect());
    }
    Ok(out)
}

/// Mirrors `return_complete_matrix(..., fast=False)`: the UCE loci shared
/// by every organism, plus per-organism loss counts (loci total minus loci
/// that organism matched at all, independent of the running intersection).
pub fn complete_matrix(
    organismal_matches: &HashMap<String, HashSet<String>>,
    organisms: &[String],
    uces: &HashSet<String>,
) -> (HashSet<String>, HashMap<String, usize>) {
    let total_loci = uces.len();
    let mut shared = uces.clone();
    let mut losses = HashMap::new();
    for organism in organisms {
        let om = organismal_matches
            .get(organism)
            .cloned()
            .unwrap_or_default();
        shared = shared.intersection(&om).cloned().collect();
        losses.insert(organism.clone(), total_loci.saturating_sub(om.len()));
    }
    (shared, losses)
}

/// Mirrors `return_incomplete_matrix`: every UCE locus matched by *any*
/// organism, intersected with the full locus set (a no-op in practice,
/// kept for parity).
pub fn incomplete_matrix(
    organismal_matches: &HashMap<String, HashSet<String>>,
    uces: &HashSet<String>,
) -> HashSet<String> {
    let all: HashSet<String> = organismal_matches.values().flatten().cloned().collect();
    all.intersection(uces).cloned().collect()
}

/// Format the `[Organisms]` / `[Loci]` output config, both lists sorted.
/// Mirrors the `main()` output-writing block's `"\n".join(sorted(...))`.
pub fn format_output(organisms: &[String], uces: &HashSet<String>) -> String {
    let mut orgs: Vec<&String> = organisms.iter().collect();
    orgs.sort();
    let mut loci: Vec<&String> = uces.iter().collect();
    loci.sort();
    format!(
        "[Organisms]\n{}\n[Loci]\n{}\n",
        orgs.iter()
            .map(|s| s.as_str())
            .collect::<Vec<_>>()
            .join("\n"),
        loci.iter()
            .map(|s| s.as_str())
            .collect::<Vec<_>>()
            .join("\n"),
    )
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn parses_allow_no_value_style_sections() {
        let dir = std::env::temp_dir().join("phyluce-assembly-conf-tests");
        std::fs::create_dir_all(&dir).unwrap();
        let path = dir.join("taxon-set.conf");
        std::fs::write(&path, "[all]\nalligator_mississippiensis\ngallus-gallus\n").unwrap();
        let cfg = read_taxon_list_config(&path).unwrap();
        assert_eq!(
            cfg.get("all").unwrap(),
            &vec![
                "alligator_mississippiensis".to_string(),
                "gallus-gallus".to_string()
            ]
        );
        let names = names_from_config(&cfg, "all").unwrap();
        assert_eq!(names, vec!["alligator_mississippiensis", "gallus_gallus"]);
    }

    #[test]
    fn excludes_are_applied() {
        let mut cfg = HashMap::new();
        cfg.insert(
            "all".to_string(),
            vec!["a".to_string(), "b".to_string(), "c".to_string()],
        );
        cfg.insert("Excludes".to_string(), vec!["b".to_string()]);
        let taxa = taxa_from_config(&cfg, "all").unwrap();
        assert_eq!(taxa, vec!["a".to_string(), "c".to_string()]);
    }

    #[test]
    fn complete_matrix_intersects_all_organisms() {
        let mut organismal = HashMap::new();
        organismal.insert(
            "a".to_string(),
            HashSet::from(["uce-1".to_string(), "uce-2".to_string()]),
        );
        organismal.insert("b".to_string(), HashSet::from(["uce-1".to_string()]));
        let uces = HashSet::from(["uce-1".to_string(), "uce-2".to_string()]);
        let organisms = vec!["a".to_string(), "b".to_string()];
        let (shared, losses) = complete_matrix(&organismal, &organisms, &uces);
        assert_eq!(shared, HashSet::from(["uce-1".to_string()]));
        assert_eq!(losses["a"], 0);
        assert_eq!(losses["b"], 1);
    }

    #[test]
    fn incomplete_matrix_unions_all_organisms() {
        let mut organismal = HashMap::new();
        organismal.insert("a".to_string(), HashSet::from(["uce-1".to_string()]));
        organismal.insert("b".to_string(), HashSet::from(["uce-2".to_string()]));
        let uces = HashSet::from([
            "uce-1".to_string(),
            "uce-2".to_string(),
            "uce-3".to_string(),
        ]);
        let result = incomplete_matrix(&organismal, &uces);
        assert_eq!(
            result,
            HashSet::from(["uce-1".to_string(), "uce-2".to_string()])
        );
    }

    #[test]
    fn formats_sorted_output() {
        let organisms = vec!["b".to_string(), "a".to_string()];
        let uces = HashSet::from(["uce-2".to_string(), "uce-1".to_string()]);
        let out = format_output(&organisms, &uces);
        assert_eq!(out, "[Organisms]\na\nb\n[Loci]\nuce-1\nuce-2\n");
    }
}
