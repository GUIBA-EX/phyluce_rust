//! phyluce-assembly: contig<->probe matching pipeline (mirrors
//! `phyluce_assembly_match_contigs_to_probes` field-for-field, including its
//! SQLite schema and CSV summary format).

use std::collections::{HashMap, HashSet};
use std::path::Path;

use regex::Regex;

use phyluce_io::lastz::LastzMatch;

pub mod explode;
pub mod get_fastas;
pub mod match_counts;
pub mod raw_reads;

#[derive(Debug, thiserror::Error)]
pub enum MatchError {
    #[error("{0}")]
    Io(#[from] std::io::Error),
    #[error("{0}")]
    Fasta(#[from] phyluce_io::FastaError),
    #[error("{0}")]
    Lastz(#[from] phyluce_io::lastz::LastzError),
    #[error("{0}")]
    Sqlite(#[from] rusqlite::Error),
    #[error("no regex match for probe name in header: {0:?}")]
    NoProbeNameMatch(String),
    #[error("no regex match for contig name in header: {0:?}")]
    NoContigNameMatch(String),
    #[error(
        "taxon name '{name}' contains or begins with an illegal character: '{illegal}'. \
         Use only letters, numbers (after a letter), and underscores"
    )]
    IllegalTaxonName { name: String, illegal: char },
}

/// Extract the locus/probe name from a probe or lastz `name2` header using
/// the (compiled) probe regex, e.g. `^(uce-\d+)(?:_p\d+.*)`. Mirrors
/// `new_get_probe_name`.
pub fn extract_probe_name(header: &str, regex: &Regex) -> Result<String, MatchError> {
    regex
        .captures(header)
        .and_then(|c| c.get(1))
        .map(|m| m.as_str().to_string())
        .ok_or_else(|| MatchError::NoProbeNameMatch(header.to_string()))
}

/// Extract the contig/node name from a lastz `name1` header using the
/// case-insensitive, anchored `[headers]` regex from phyluce.conf. Mirrors
/// `get_contig_name`.
pub fn extract_contig_name(header: &str, header_regex: &Regex) -> Result<String, MatchError> {
    header_regex
        .captures(header)
        .and_then(|c| c.get(1))
        .map(|m| m.as_str().to_string())
        .ok_or_else(|| MatchError::NoContigNameMatch(header.to_string()))
}

/// Build the case-insensitive, anchored contig-header regex from the
/// `|`-joined `[headers]` fragments, matching
/// `re.search("^({}).*".format(contig_header_string), header, flags=re.I)`.
pub fn contig_header_regex(joined_header_fragments: &str) -> Result<Regex, regex::Error> {
    Regex::new(&format!("(?i)^({}).*", joined_header_fragments))
}

/// Count FASTA records in a file by counting header lines -- mirrors
/// `contig_count` (a raw `>`-prefix line count, not a full FASTA parse).
pub fn contig_count(path: &Path) -> std::io::Result<usize> {
    use std::io::BufRead;
    let f = std::fs::File::open(path)?;
    let reader = std::io::BufReader::new(f);
    Ok(reader
        .lines()
        .map_while(Result::ok)
        .filter(|l| l.starts_with('>'))
        .count())
}

const ILLEGAL_TAXON_CHARS: &str = ".+:\"'-?!*@%^&#=/\\";

/// Parse organism/taxon names from contig FASTA file paths (basename before
/// the first `.`, `-` replaced with `_`), then validate that no name begins
/// with a digit or contains a character SQLite can't use as an identifier.
/// Mirrors `get_organism_names_from_fasta_files`.
pub fn organism_names_from_fasta_paths(
    paths: &[std::path::PathBuf],
) -> Result<Vec<String>, Vec<MatchError>> {
    let names: Vec<String> = paths
        .iter()
        .map(|p| {
            let stem = p
                .file_name()
                .and_then(|n| n.to_str())
                .unwrap_or("")
                .split('.')
                .next()
                .unwrap_or("")
                .to_string();
            stem.replace('-', "_")
        })
        .collect();

    let mut errors = Vec::new();
    for name in &names {
        if let Some(c) = name.chars().next() {
            if c.is_ascii_digit() {
                errors.push(MatchError::IllegalTaxonName {
                    name: name.clone(),
                    illegal: c,
                });
                continue;
            }
        }
        if let Some(c) = name.chars().find(|c| ILLEGAL_TAXON_CHARS.contains(*c)) {
            errors.push(MatchError::IllegalTaxonName {
                name: name.clone(),
                illegal: c,
            });
        }
    }
    if errors.is_empty() {
        Ok(names)
    } else {
        Err(errors)
    }
}

/// Result of processing one taxon's LASTZ alignment against the probe set,
/// before duplicate contigs/loci are dropped.
#[derive(Debug, Default)]
pub struct TaxonMatches {
    /// contig name -> set of UCE loci it hit
    pub matches: HashMap<String, HashSet<String>>,
    /// UCE locus -> set of strands it was hit on (usually a single strand)
    pub orientation: HashMap<String, HashSet<String>>,
    /// UCE locus -> set of contigs that hit it
    pub revmatches: HashMap<String, HashSet<String>>,
    /// UCE loci excluded because they're flagged as probe duplicates
    pub probe_dupes: HashSet<String>,
}

/// Mirrors the `main()` per-taxon LASTZ-result loop: classify each match by
/// contig/locus, tracking probe duplicates separately.
pub fn process_taxon_lastz(
    lastz_matches: &[LastzMatch],
    probe_regex: &Regex,
    header_regex: &Regex,
    dupes: &HashSet<String>,
    dupefile_active: bool,
) -> Result<TaxonMatches, MatchError> {
    let mut result = TaxonMatches::default();
    for m in lastz_matches {
        let contig_name = extract_contig_name(&m.name1, header_regex)?;
        let uce_name = extract_probe_name(&m.name2, probe_regex)?;
        if dupefile_active && dupes.contains(&uce_name) {
            result.probe_dupes.insert(uce_name);
        } else {
            result
                .matches
                .entry(contig_name.clone())
                .or_default()
                .insert(uce_name.clone());
            result
                .orientation
                .entry(uce_name.clone())
                .or_default()
                .insert(m.strand2.clone());
            result
                .revmatches
                .entry(uce_name)
                .or_default()
                .insert(contig_name);
        }
    }
    Ok(result)
}

/// Contigs that matched more than one distinct UCE locus. Mirrors
/// `check_contigs_for_dupes`.
pub fn contigs_matching_multiple_uces(
    matches: &HashMap<String, HashSet<String>>,
) -> HashSet<String> {
    matches
        .iter()
        .filter(|(_, uces)| uces.len() > 1)
        .map(|(node, _)| node.clone())
        .collect()
}

/// (contigs, UCE loci) where a UCE locus matched more than one contig.
/// Mirrors `check_loci_for_dupes`.
pub fn loci_matching_multiple_contigs(
    revmatches: &HashMap<String, HashSet<String>>,
) -> (HashSet<String>, HashSet<String>) {
    let mut dupe_contigs = HashSet::new();
    let mut dupe_uces = HashSet::new();
    for (uce, contigs) in revmatches {
        if contigs.len() > 1 {
            dupe_contigs.extend(contigs.iter().cloned());
            dupe_uces.insert(uce.clone());
        }
    }
    (dupe_contigs, dupe_uces)
}

/// Self-to-self LASTZ duplicate probe detection, keyed by the
/// regex-extracted probe/locus name (not the raw header) -- mirrors the
/// script-local `get_dupes` (distinct from `helpers.get_dupes`, which keys
/// on a split character instead of a regex).
pub fn get_probe_dupes(
    lastz_matches: &[LastzMatch],
    probe_regex: &Regex,
) -> Result<HashSet<String>, MatchError> {
    let mut matches: HashMap<String, Vec<String>> = HashMap::new();
    for m in lastz_matches {
        let target = extract_probe_name(&m.name1, probe_regex)?;
        let query = extract_probe_name(&m.name2, probe_regex)?;
        matches.entry(target).or_default().push(query);
    }
    let mut dupes = HashSet::new();
    for (k, v) in &matches {
        if v.len() > 1 {
            for i in v {
                if i != k {
                    dupes.insert(k.clone());
                    dupes.insert(i.clone());
                }
            }
        } else if &v[0] != k {
            dupes.insert(k.clone());
        }
    }
    Ok(dupes.into_iter().map(|d| d.to_lowercase()).collect())
}

pub mod db {
    use super::*;
    use phyluce_io::sql::ident;
    use rusqlite::Connection;

    /// Create `probe.matches.sqlite`'s `matches`/`match_map` tables (one
    /// text column per organism, `uce` as the primary key) and seed every
    /// UCE locus as a row. Mirrors `create_probe_database`.
    pub fn create_probe_database(
        path: &Path,
        organisms: &[String],
        uces: &[String],
    ) -> rusqlite::Result<Connection> {
        let conn = Connection::open(path)?;
        conn.execute_batch("PRAGMA foreign_keys = ON;")?;
        let columns = organisms
            .iter()
            .map(|o| format!("{} text", ident(o)))
            .collect::<Vec<_>>()
            .join(",");
        conn.execute(
            &format!("CREATE TABLE matches (uce text primary key, {columns})"),
            [],
        )?;
        conn.execute(
            &format!("CREATE TABLE match_map (uce text primary key, {columns})"),
            [],
        )?;
        let tx = conn.unchecked_transaction()?;
        {
            let mut ins_matches = tx.prepare("INSERT INTO matches(uce) VALUES (?1)")?;
            let mut ins_map = tx.prepare("INSERT INTO match_map(uce) VALUES (?1)")?;
            for uce in uces {
                ins_matches.execute([uce])?;
                ins_map.execute([uce])?;
            }
        }
        tx.commit()?;
        Ok(conn)
    }

    /// Mirrors `store_lastz_results_in_db`: for every (contig -> single UCE)
    /// pair left after duplicate filtering, mark the match and record
    /// `contig(strand)` in `match_map`.
    pub fn store_lastz_results(
        conn: &Connection,
        matches: &HashMap<String, HashSet<String>>,
        orientation: &HashMap<String, HashSet<String>>,
        critter: &str,
    ) -> rusqlite::Result<()> {
        for (contig, uces) in matches {
            assert_eq!(uces.len(), 1, "More than one match");
            let uce = uces.iter().next().unwrap();
            let column = ident(critter);
            conn.execute(
                &format!("UPDATE matches SET {column} = 1 WHERE uce = ?1"),
                [uce],
            )?;
            let strand = orientation
                .get(uce)
                .and_then(|s| s.iter().next())
                .cloned()
                .unwrap_or_default();
            let orient_value = format!("{contig}({strand})");
            conn.execute(
                &format!("UPDATE match_map SET {column} = ?1 WHERE uce = ?2"),
                rusqlite::params![orient_value, uce],
            )?;
        }
        Ok(())
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn probe_regex() -> Regex {
        Regex::new(r"^(uce-\d+)(?:_p\d+.*)").unwrap()
    }

    #[test]
    fn extracts_probe_name() {
        let re = probe_regex();
        assert_eq!(
            extract_probe_name("uce-553_p1 |source:x", &re).unwrap(),
            "uce-553"
        );
    }

    #[test]
    fn extracts_contig_name_case_insensitively() {
        let re = contig_header_regex(r"node_\d+").unwrap();
        assert_eq!(
            extract_contig_name("NODE_1_length_500_cov_10", &re).unwrap(),
            "NODE_1"
        );
    }

    #[test]
    fn organism_names_replace_dash_with_underscore() {
        let paths = vec![std::path::PathBuf::from(
            "/x/alligator-mississippiensis.contigs.fasta",
        )];
        let names = organism_names_from_fasta_paths(&paths).unwrap();
        assert_eq!(names, vec!["alligator_mississippiensis"]);
    }

    #[test]
    fn organism_names_reject_leading_digit() {
        let paths = vec![std::path::PathBuf::from("/x/1bad.contigs.fasta")];
        let err = organism_names_from_fasta_paths(&paths).unwrap_err();
        assert_eq!(err.len(), 1);
    }

    #[test]
    fn detects_contigs_matching_multiple_uces() {
        let mut matches = HashMap::new();
        matches.insert(
            "NODE_1".to_string(),
            HashSet::from(["uce-1".to_string(), "uce-2".to_string()]),
        );
        matches.insert("NODE_2".to_string(), HashSet::from(["uce-3".to_string()]));
        let dupes = contigs_matching_multiple_uces(&matches);
        assert_eq!(dupes, HashSet::from(["NODE_1".to_string()]));
    }

    #[test]
    fn detects_loci_matching_multiple_contigs() {
        let mut revmatches = HashMap::new();
        revmatches.insert(
            "uce-1".to_string(),
            HashSet::from(["NODE_1".to_string(), "NODE_2".to_string()]),
        );
        revmatches.insert("uce-2".to_string(), HashSet::from(["NODE_3".to_string()]));
        let (dupe_contigs, dupe_uces) = loci_matching_multiple_contigs(&revmatches);
        assert_eq!(
            dupe_contigs,
            HashSet::from(["NODE_1".to_string(), "NODE_2".to_string()])
        );
        assert_eq!(dupe_uces, HashSet::from(["uce-1".to_string()]));
    }

    #[test]
    fn db_roundtrip_creates_schema_and_seeds_rows() {
        let dir = std::env::temp_dir().join("phyluce-assembly-db-tests");
        std::fs::create_dir_all(&dir).unwrap();
        let db_path = dir.join(format!("test-{}.sqlite", std::process::id()));
        let _ = std::fs::remove_file(&db_path);

        let organisms = vec!["taxon_a".to_string(), "taxon_b".to_string()];
        let uces = vec!["uce-1".to_string(), "uce-2".to_string()];
        let conn = db::create_probe_database(&db_path, &organisms, &uces).unwrap();

        let mut matches = HashMap::new();
        matches.insert("NODE_1".to_string(), HashSet::from(["uce-1".to_string()]));
        let mut orientation = HashMap::new();
        orientation.insert("uce-1".to_string(), HashSet::from(["+".to_string()]));
        db::store_lastz_results(&conn, &matches, &orientation, "taxon_a").unwrap();

        let value: Option<String> = conn
            .query_row(
                "SELECT taxon_a FROM match_map WHERE uce = 'uce-1'",
                [],
                |r| r.get(0),
            )
            .unwrap();
        assert_eq!(value.as_deref(), Some("NODE_1(+)"));

        let flag: Option<String> = conn
            .query_row("SELECT taxon_a FROM matches WHERE uce = 'uce-1'", [], |r| {
                r.get(0)
            })
            .unwrap();
        assert_eq!(flag.as_deref(), Some("1"));

        drop(conn);
        let _ = std::fs::remove_file(&db_path);
    }
}
