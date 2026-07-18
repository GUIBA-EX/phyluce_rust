//! Contig-count / matrix-membership pipeline mirroring
//! `phyluce_assembly_get_match_counts`, including complete/incomplete matrix
//! generation and complete-matrix taxon-group optimization.

use std::collections::{BTreeMap, HashMap, HashSet};
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
    #[error("optimization group size must be between 1 and {taxa}, got {size}")]
    InvalidOptimizationGroupSize { size: usize, taxa: usize },
    #[error("--samples must be greater than zero")]
    InvalidSampleCount,
    #[error(
        "--sample-size must be between 1 and one less than the number of taxa ({taxa}), got {size}"
    )]
    InvalidSampleSize { size: usize, taxa: usize },
    #[error("taxon {taxon:?} is not a column in {table}")]
    UnknownTaxonColumn { taxon: String, table: &'static str },
}

#[derive(Clone, Debug, Eq, PartialEq)]
pub struct OptimizedGroup {
    pub organisms: Vec<String>,
    pub uces: HashSet<String>,
}

impl OptimizedGroup {
    pub fn locus_count(&self) -> usize {
        self.uces.len()
    }
}

#[derive(Clone, Debug, Eq, PartialEq)]
pub struct SampleOptimization {
    pub best: OptimizedGroup,
    pub counts: Vec<(usize, usize)>,
    pub missing_counts: BTreeMap<String, usize>,
}

/// A `[section]`-per-line conf file, `allow_no_value=True`-style: each
/// non-empty, non-comment line under a section is an item name (an
/// optional `:value`/`=value` suffix is ignored). Mirrors
/// `configparser.RawConfigParser(allow_no_value=True)` as used for
/// `--taxon-list-config`.
pub fn read_taxon_list_config(path: &Path) -> std::io::Result<HashMap<String, Vec<String>>> {
    let text = std::fs::read_to_string(path)?;
    let ini = phyluce_config::Ini::parse_allow_no_value(&text);
    Ok(ini
        .section_names()
        .map(|section| {
            let entries = ini
                .entries(section)
                .unwrap_or_default()
                .iter()
                .map(|(key, _)| key.clone())
                .collect();
            (section.to_string(), entries)
        })
        .collect())
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
) -> Result<HashMap<String, HashSet<String>>, MatchCountError> {
    let mut out = HashMap::new();
    for organism in organisms {
        let (table, column) = if let Some(stripped) = organism.strip_suffix('*') {
            ("extended.matches", stripped.to_string())
        } else {
            ("matches", organism.clone())
        };
        ensure_taxon_column(conn, table, &column, organism)?;
        let map_table = if organism.ends_with('*') {
            "extended.match_map"
        } else {
            "match_map"
        };
        ensure_taxon_column(conn, map_table, &column, organism)?;
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

fn ensure_taxon_column(
    conn: &Connection,
    table: &'static str,
    column: &str,
    taxon: &str,
) -> Result<(), MatchCountError> {
    let pragma = match table {
        "matches" => "PRAGMA table_info(matches)",
        "match_map" => "PRAGMA table_info(match_map)",
        "extended.matches" => "PRAGMA extended.table_info(matches)",
        "extended.match_map" => "PRAGMA extended.table_info(match_map)",
        _ => unreachable!("only match tables are validated"),
    };
    let mut statement = conn.prepare(pragma)?;
    let columns: HashSet<String> = statement
        .query_map([], |row| row.get(1))?
        .collect::<Result<_, _>>()?;
    if columns.contains(column) {
        Ok(())
    } else {
        Err(MatchCountError::UnknownTaxonColumn {
            taxon: taxon.to_string(),
            table,
        })
    }
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

fn intersect_matches(
    shared: &HashSet<String>,
    matches: Option<&HashSet<String>>,
) -> HashSet<String> {
    let Some(matches) = matches else {
        return HashSet::new();
    };
    if shared.len() <= matches.len() {
        shared
            .iter()
            .filter(|uce| matches.contains(*uce))
            .cloned()
            .collect()
    } else {
        matches
            .iter()
            .filter(|uce| shared.contains(*uce))
            .cloned()
            .collect()
    }
}

fn search_combinations(
    organismal_matches: &HashMap<String, HashSet<String>>,
    organisms: &[String],
    start: usize,
    remaining: usize,
    chosen: &mut Vec<usize>,
    shared: &HashSet<String>,
    best: &mut Option<OptimizedGroup>,
) {
    if best
        .as_ref()
        .is_some_and(|current| shared.len() <= current.locus_count())
    {
        return;
    }
    if remaining == 0 {
        *best = Some(OptimizedGroup {
            organisms: chosen
                .iter()
                .map(|&index| organisms[index].clone())
                .collect(),
            uces: shared.clone(),
        });
        return;
    }

    let last_start = organisms.len() - remaining;
    for index in start..=last_start {
        let next_shared = intersect_matches(shared, organismal_matches.get(&organisms[index]));
        chosen.push(index);
        search_combinations(
            organismal_matches,
            organisms,
            index + 1,
            remaining - 1,
            chosen,
            &next_shared,
            best,
        );
        chosen.pop();
    }
}

/// Return the first taxon combination of `size` with the largest complete
/// matrix. Strict `>` replacement preserves the legacy first-combination tie
/// rule.
pub fn optimize_group_for_size(
    organismal_matches: &HashMap<String, HashSet<String>>,
    organisms: &[String],
    uces: &HashSet<String>,
    size: usize,
) -> Result<OptimizedGroup, MatchCountError> {
    if size == 0 || size > organisms.len() {
        return Err(MatchCountError::InvalidOptimizationGroupSize {
            size,
            taxa: organisms.len(),
        });
    }
    let mut best = None;
    search_combinations(
        organismal_matches,
        organisms,
        0,
        size,
        &mut Vec::with_capacity(size),
        uces,
        &mut best,
    );
    Ok(best.expect("a valid group size always has at least one combination"))
}

struct SplitMix64 {
    state: u64,
}

impl SplitMix64 {
    fn new(seed: u64) -> Self {
        Self { state: seed }
    }

    fn next_u64(&mut self) -> u64 {
        self.state = self.state.wrapping_add(0x9E3779B97F4A7C15);
        let mut value = self.state;
        value = (value ^ (value >> 30)).wrapping_mul(0xBF58476D1CE4E5B9);
        value = (value ^ (value >> 27)).wrapping_mul(0x94D049BB133111EB);
        value ^ (value >> 31)
    }

    fn index(&mut self, upper: usize) -> usize {
        debug_assert!(upper > 0);
        let upper = upper as u64;
        let threshold = upper.wrapping_neg() % upper;
        loop {
            let value = self.next_u64();
            if value >= threshold {
                return (value % upper) as usize;
            }
        }
    }
}

fn sample_without_replacement(
    organisms: &[String],
    count: usize,
    rng: &mut SplitMix64,
) -> Vec<String> {
    let mut pool = organisms.to_vec();
    for index in 0..count {
        let selected = index + rng.index(pool.len() - index);
        pool.swap(index, selected);
    }
    pool.truncate(count);
    pool
}

/// Mirror the legacy sampling objective: draw `sample_size + 1` taxa, then
/// retain the `sample_size`-taxon subset with the largest complete matrix.
pub fn sample_optimized_groups(
    organismal_matches: &HashMap<String, HashSet<String>>,
    organisms: &[String],
    uces: &HashSet<String>,
    samples: usize,
    sample_size: usize,
    seed: u64,
) -> Result<SampleOptimization, MatchCountError> {
    if samples == 0 {
        return Err(MatchCountError::InvalidSampleCount);
    }
    if sample_size == 0 || sample_size >= organisms.len() {
        return Err(MatchCountError::InvalidSampleSize {
            size: sample_size,
            taxa: organisms.len(),
        });
    }

    let mut rng = SplitMix64::new(seed);
    let mut best: Option<OptimizedGroup> = None;
    let mut counts = Vec::with_capacity(samples);
    let mut missing_counts = BTreeMap::new();

    for _ in 0..samples {
        let sampled = sample_without_replacement(organisms, sample_size + 1, &mut rng);
        let group = optimize_group_for_size(organismal_matches, &sampled, uces, sample_size)?;
        counts.push((sample_size, group.locus_count()));

        let selected: HashSet<&str> = group.organisms.iter().map(String::as_str).collect();
        for organism in organisms {
            if !selected.contains(organism.as_str()) {
                *missing_counts.entry(organism.clone()).or_insert(0) += 1;
            }
        }

        if best
            .as_ref()
            .is_none_or(|current| group.locus_count() > current.locus_count())
        {
            best = Some(group);
        }
    }

    Ok(SampleOptimization {
        best: best.expect("positive sample count always produces a group"),
        counts,
        missing_counts,
    })
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
        std::fs::write(
            &path,
            "[all]\n; an ordinary INI comment\nalligator_mississippiensis\ngallus-gallus\n",
        )
        .unwrap();
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
    fn rejects_unknown_taxon_columns() {
        let conn = Connection::open_in_memory().unwrap();
        conn.execute_batch(
            "CREATE TABLE matches (uce TEXT PRIMARY KEY, known INTEGER);
             CREATE TABLE match_map (uce TEXT PRIMARY KEY, known TEXT);",
        )
        .unwrap();
        let error = matches_by_organism(&conn, &["missing".to_string()]).unwrap_err();
        assert!(matches!(
            error,
            MatchCountError::UnknownTaxonColumn { ref taxon, table: "matches" } if taxon == "missing"
        ));
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

    fn optimization_fixture() -> (
        HashMap<String, HashSet<String>>,
        Vec<String>,
        HashSet<String>,
    ) {
        let organisms = ["a", "b", "c", "d"].map(str::to_string).to_vec();
        let organismal = HashMap::from([
            (
                "a".to_string(),
                ["uce-1", "uce-2", "uce-3", "uce-4"]
                    .map(str::to_string)
                    .into_iter()
                    .collect(),
            ),
            (
                "b".to_string(),
                ["uce-1", "uce-2", "uce-3"]
                    .map(str::to_string)
                    .into_iter()
                    .collect(),
            ),
            (
                "c".to_string(),
                ["uce-1", "uce-2", "uce-4"]
                    .map(str::to_string)
                    .into_iter()
                    .collect(),
            ),
            (
                "d".to_string(),
                ["uce-1", "uce-4"].map(str::to_string).into_iter().collect(),
            ),
        ]);
        let uces = ["uce-1", "uce-2", "uce-3", "uce-4"]
            .map(str::to_string)
            .into_iter()
            .collect();
        (organismal, organisms, uces)
    }

    #[test]
    fn optimization_preserves_first_combination_on_ties() {
        let (organismal, organisms, uces) = optimization_fixture();
        let best = optimize_group_for_size(&organismal, &organisms, &uces, 2).unwrap();
        assert_eq!(best.organisms, vec!["a", "b"]);
        assert_eq!(best.locus_count(), 3);

        let best = optimize_group_for_size(&organismal, &organisms, &uces, 3).unwrap();
        assert_eq!(best.organisms, vec!["a", "b", "c"]);
        assert_eq!(best.locus_count(), 2);
    }

    #[test]
    fn sampled_optimization_is_reproducible_and_counts_missing_taxa() {
        let (organismal, organisms, uces) = optimization_fixture();
        let first = sample_optimized_groups(&organismal, &organisms, &uces, 20, 2, 42).unwrap();
        let second = sample_optimized_groups(&organismal, &organisms, &uces, 20, 2, 42).unwrap();
        assert_eq!(first, second);
        assert_eq!(first.counts.len(), 20);
        assert_eq!(first.best.organisms.len(), 2);
        assert_eq!(first.missing_counts.values().sum::<usize>(), 40);
    }

    #[test]
    fn sampled_optimization_rejects_invalid_ranges() {
        let (organismal, organisms, uces) = optimization_fixture();
        assert!(sample_optimized_groups(&organismal, &organisms, &uces, 0, 2, 1).is_err());
        assert!(sample_optimized_groups(&organismal, &organisms, &uces, 1, 0, 1).is_err());
        assert!(sample_optimized_groups(&organismal, &organisms, &uces, 1, 4, 1).is_err());
        assert!(optimize_group_for_size(&organismal, &organisms, &uces, 0).is_err());
    }

    #[test]
    fn formats_sorted_output() {
        let organisms = vec!["b".to_string(), "a".to_string()];
        let uces = HashSet::from(["uce-2".to_string(), "uce-1".to_string()]);
        let out = format_output(&organisms, &uces);
        assert_eq!(out, "[Organisms]\na\nb\n[Loci]\nuce-1\nuce-2\n");
    }
}
