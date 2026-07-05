//! Locus-extraction pipeline mirroring
//! `phyluce_assembly_get_fastas_from_match_counts`: for each taxon in a
//! match-count-output config, find its contig FASTA, pull out every
//! matched UCE locus, rename/reorient/clean it, and write a monolithic
//! output FASTA.

use std::collections::HashMap;
use std::path::{Path, PathBuf};

use phyluce_io::sql::{ident, qualified_ident};
use regex::Regex;
use rusqlite::Connection;

#[derive(Debug, thiserror::Error)]
pub enum GetFastasError {
    #[error("{0}")]
    Io(#[from] std::io::Error),
    #[error("{0}")]
    Sqlite(#[from] rusqlite::Error),
    #[error("could not parse contig name/strand from match_map value: {0:?}")]
    BadNodeValue(String),
    #[error("complete matrices should have no missing data (organism {organism:?}, uce {uce:?})")]
    UnexpectedMissingData { organism: String, uce: String },
    #[error("cannot find a fasta file for {name} with any of the extensions ({extensions})")]
    ContigFileNotFound { name: String, extensions: String },
}

/// One resolved UCE match: the UCE locus name and the strand the probe hit
/// on (used to decide whether to reverse-complement the extracted contig).
#[derive(Debug, Clone)]
pub struct NodeMatch {
    pub uce: String,
    pub strand: char,
}

const CONTIG_EXTENSIONS: &[&str] = &[
    ".fa",
    ".fasta",
    ".contigs.fasta",
    ".contigs.fa",
    ".gz",
    ".fasta.gz",
    ".fa.gz",
];

/// Mirrors `get_nodes_for_uces`: for the given organism, resolve every
/// requested UCE locus to the contig node (lowercased) and strand that
/// matched it, splitting missing/`NULL` matches into `missing` when
/// `notstrict` (an incomplete matrix); otherwise a `NULL` match is an
/// error (a complete matrix should never have one).
pub fn get_nodes_for_uces(
    conn: &Connection,
    organism: &str,
    uces: &[String],
    extend: bool,
    notstrict: bool,
    header_with_strand_regex: &Regex,
) -> Result<(HashMap<String, NodeMatch>, Vec<String>), GetFastasError> {
    if uces.is_empty() {
        return Ok((HashMap::new(), Vec::new()));
    }
    let table = if extend {
        "extended.match_map"
    } else {
        "match_map"
    };
    let in_list = vec!["?"; uces.len()].join(",");
    let column = ident(organism);
    let table = qualified_ident(table);
    let query = format!("SELECT lower({column}), uce FROM {table} where uce in ({in_list})");
    let mut stmt = conn.prepare(&query)?;
    let rows: Vec<(Option<String>, String)> = stmt
        .query_map(rusqlite::params_from_iter(uces.iter()), |r| {
            Ok((r.get::<_, Option<String>>(0)?, r.get::<_, String>(1)?))
        })?
        .collect::<Result<_, _>>()?;

    let mut node_dict = HashMap::new();
    let mut missing = Vec::new();
    for (node_value, uce) in rows {
        match node_value {
            Some(v) => {
                let caps = header_with_strand_regex
                    .captures(&v)
                    .ok_or_else(|| GetFastasError::BadNodeValue(v.clone()))?;
                let node = caps.get(1).unwrap().as_str().to_string();
                let strand = caps.get(2).unwrap().as_str().chars().next().unwrap();
                node_dict.insert(node, NodeMatch { uce, strand });
            }
            None if notstrict => missing.push(uce),
            None => {
                return Err(GetFastasError::UnexpectedMissingData {
                    organism: organism.to_string(),
                    uce,
                })
            }
        }
    }
    Ok((node_dict, missing))
}

/// Build the `^({header_fragments})\(([+-])\)` (case-insensitive) regex
/// used to split a lowercased `match_map` value like `node_1_length_500(+)`
/// into (contig node, strand). Mirrors `get_nodes_for_uces`'s inline regex.
pub fn node_with_strand_regex(header_fragments: &str) -> Result<Regex, regex::Error> {
    Regex::new(&format!(r"(?i)^({header_fragments})\(([+-])\)"))
}

/// Mirrors `find_file`: try `name` and `name` with `-`/`_` swapped, across
/// a fixed list of extensions, also trying an all-lowercase path.
pub fn find_contig_file(contigs_dir: &Path, name: &str) -> Result<PathBuf, GetFastasError> {
    let alt_name = if name.contains('-') {
        name.replace('-', "_")
    } else {
        name.replace('_', "-")
    };
    for ext in CONTIG_EXTENSIONS {
        for candidate_name in [name, alt_name.as_str()] {
            let candidate = contigs_dir.join(format!("{candidate_name}{ext}"));
            if candidate.is_file() {
                return Ok(candidate);
            }
            let lower = contigs_dir.join(format!("{candidate_name}{ext}").to_lowercase());
            if lower.is_file() {
                return Ok(lower);
            }
        }
    }
    Err(GetFastasError::ContigFileNotFound {
        name: name.to_string(),
        extensions: CONTIG_EXTENSIONS.join(", "),
    })
}

/// Runs of these characters (regardless of length) are deleted from the
/// raw sequence -- mirrors `re.compile("[N,n]{1,21}")` applied globally
/// (chunked non-overlapping `{1,21}` matches still consume an arbitrarily
/// long run in full, so an unbounded `+` is equivalent here).
fn strip_n_runs(seq: &str) -> (String, bool) {
    let matched = seq.contains(['N', 'n', ',']);
    let mut out = String::with_capacity(seq.len());
    for c in seq.chars() {
        if c != 'N' && c != 'n' && c != ',' {
            out.push(c);
        }
    }
    (out, matched)
}

/// Mirrors `replace_and_remove_bases`: delete embedded N/n/`,` runs, then
/// strip leading/trailing runs of lowercase `acgtn` (soft-masked,
/// low-coverage bases from some assemblers).
pub fn clean_sequence(seq: &str) -> (String, bool) {
    let (deambiguated, replaced) = strip_n_runs(seq);
    let trimmed = deambiguated
        .trim_start_matches(|c: char| "acgtn".contains(c))
        .trim_end_matches(|c: char| "acgtn".contains(c));
    (trimmed.to_string(), replaced)
}

/// IUPAC DNA complement (upper- and lower-case), matching Biopython's
/// `ambiguous_dna_complement` table; unrecognized characters (gaps, etc.)
/// pass through unchanged.
fn complement_base(c: char) -> char {
    match c {
        'A' => 'T',
        'T' => 'A',
        'C' => 'G',
        'G' => 'C',
        'M' => 'K',
        'K' => 'M',
        'R' => 'Y',
        'Y' => 'R',
        'W' => 'W',
        'S' => 'S',
        'V' => 'B',
        'B' => 'V',
        'H' => 'D',
        'D' => 'H',
        'X' => 'X',
        'N' => 'N',
        'a' => 't',
        't' => 'a',
        'c' => 'g',
        'g' => 'c',
        'm' => 'k',
        'k' => 'm',
        'r' => 'y',
        'y' => 'r',
        'w' => 'w',
        's' => 's',
        'v' => 'b',
        'b' => 'v',
        'h' => 'd',
        'd' => 'h',
        'x' => 'x',
        'n' => 'n',
        other => other,
    }
}

/// Mirrors `seq.seq.reverse_complement()`.
pub fn reverse_complement(seq: &str) -> String {
    seq.chars().rev().map(complement_base).collect()
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn strips_internal_and_edge_n_runs() {
        let (cleaned, replaced) = clean_sequence(
            "NNNNacgtACGTnnnNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNacgt",
        );
        assert!(replaced);
        assert_eq!(cleaned, "ACGT");
    }

    #[test]
    fn strips_only_lowercase_edges() {
        let (cleaned, replaced) = clean_sequence("acgtACGTacgt");
        assert!(!replaced);
        assert_eq!(cleaned, "ACGT");
    }

    #[test]
    fn reverse_complement_handles_iupac_ambiguity() {
        assert_eq!(reverse_complement("ACGTMRWSN"), "NSWYKACGT");
    }

    #[test]
    fn find_contig_file_tries_hyphen_and_underscore() {
        let dir = std::env::temp_dir().join("phyluce-assembly-findfile-tests");
        std::fs::create_dir_all(&dir).unwrap();
        let path = dir.join("alligator_mississippiensis.contigs.fasta");
        std::fs::write(&path, ">a\nACGT\n").unwrap();
        let found = find_contig_file(&dir, "alligator-mississippiensis").unwrap();
        assert_eq!(found, path);
    }
}
