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

/// `HashMap`/`HashSet` keyed on `ahash` instead of the standard library's
/// SipHash. Benchmarked in `tests::bench_process_taxon_lastz_end_to_end`:
/// the contig/UCE matching pipeline is dominated by regex extraction
/// (~57% of wall time), not hashing (~10%), so this buys a modest,
/// low-risk win on the hottest maps/sets rather than a dramatic one --
/// see the benchmark module doc comment for the full breakdown. Safe to
/// swap because every consumer of these collections already sorts before
/// producing output (SQLite writes are keyed `UPDATE`s, not appends; CSV/
/// report output goes through `Vec` + `.sort()`), so none of it depends on
/// iteration order.
pub type FastMap<K, V> = HashMap<K, V, ahash::RandomState>;
pub type FastSet<T> = HashSet<T, ahash::RandomState>;

/// Hand-rolled equivalents of phyluce's *default* `--regex` (probe/locus
/// name extraction) and `[headers]` (contig name extraction) patterns,
/// used as a fast path in `extract_probe_name`/`extract_contig_name` when
/// the caller's compiled regex is byte-identical to those defaults.
///
/// Both patterns are runtime-configurable (`--regex` on the CLI, `[headers]`
/// in `phyluce.conf`), so this can't just replace the regex engine --
/// anything other than an exact match on the default source falls through
/// to `Regex::captures` unchanged. Benchmarked in
/// `tests::bench_regex_vs_hand_rolled_default_pattern_ceiling`: ~12x faster
/// than `regex` on the probe pattern alone, which is worth the maintenance
/// cost *only* because every function here is a pure "does this match, and
/// if so how long is the match" check with no independent error path --
/// a `None` always defers to the real regex rather than reporting "no
/// match" on its own, so a bug in this module can only cost performance
/// (falling back to the slow, definitely-correct path), never produce a
/// wrong contig/probe name. See `tests::fast_extract_matches_regex_*` for
/// the differential fuzz tests that back that claim.
mod fast_extract {
    /// Mirrors `EasyLastz`'s / `match-contigs-to-probes`'s default
    /// `--regex` value in `phyluce-cli/src/main.rs`. Kept in sync by
    /// `tests::default_probe_regex_source_matches_cli_default`.
    pub const DEFAULT_PROBE_REGEX_SRC: &str = r"^(uce-\d+)(?:_p\d+.*)";

    /// The default `header_regex` source, derived at runtime (not
    /// hardcoded) from the packaged `config/phyluce.conf`'s `[headers]`
    /// section via `contig_header_regex`, so this never drifts from the
    /// actual shipped defaults.
    pub fn default_header_regex_source() -> &'static str {
        static SRC: std::sync::OnceLock<String> = std::sync::OnceLock::new();
        SRC.get_or_init(|| {
            phyluce_config::PhyluceConfig::load_package_only()
                .ok()
                .and_then(|cfg| cfg.get_contig_header_string())
                .and_then(|frags| super::contig_header_regex(&frags).ok())
                .map(|re| re.as_str().to_string())
                // If the packaged config can't be loaded/parsed for any
                // reason, this never equals a real `Regex::as_str()`
                // output, so the fast path just never activates.
                .unwrap_or_default()
        })
    }

    fn eat_ci_literal(b: &[u8], pos: usize, lit: &str) -> Option<usize> {
        let lit = lit.as_bytes();
        if pos + lit.len() > b.len() {
            return None;
        }
        for i in 0..lit.len() {
            if !b[pos + i].eq_ignore_ascii_case(&lit[i]) {
                return None;
            }
        }
        Some(pos + lit.len())
    }

    /// Case-*sensitive* literal match. `DEFAULT_PROBE_REGEX_SRC` has no
    /// `(?i)` flag (unlike the header pattern, which does), so `probe_name`
    /// must not fold case -- an earlier version of this used
    /// `eat_ci_literal` here and the differential fuzz test
    /// (`fast_extract_probe_name_matches_regex_oracle`) caught it
    /// producing false positives on mixed-case input like `"uCE-1_P2"`.
    fn eat_literal(b: &[u8], pos: usize, lit: &str) -> Option<usize> {
        let lit = lit.as_bytes();
        if pos + lit.len() > b.len() {
            return None;
        }
        if &b[pos..pos + lit.len()] == lit {
            Some(pos + lit.len())
        } else {
            None
        }
    }

    fn eat_digits(b: &[u8], pos: usize) -> Option<usize> {
        let mut p = pos;
        while p < b.len() && b[p].is_ascii_digit() {
            p += 1;
        }
        if p > pos {
            Some(p)
        } else {
            None
        }
    }

    /// Regex `.`: matches exactly one arbitrary byte (headers are
    /// single-line ASCII, so this never needs to special-case `\n`).
    fn eat_any_byte(b: &[u8], pos: usize) -> Option<usize> {
        if pos < b.len() {
            Some(pos + 1)
        } else {
            None
        }
    }

    /// `uce-\d+` (the captured group) followed by the required, uncaptured
    /// `_p\d+.*`. Returns the end offset of the *captured* group only.
    pub fn probe_name(header: &str) -> Option<&str> {
        let b = header.as_bytes();
        let p = eat_literal(b, 0, "uce-")?;
        let end = eat_digits(b, p)?;
        let p = eat_literal(b, end, "_p")?;
        eat_digits(b, p)?;
        Some(&header[..end])
    }

    // One function per `|`-separated alternative in the default
    // `[headers]` fragments, in the exact left-to-right order they appear
    // in `config/phyluce.conf` -- regex alternation is leftmost-first, not
    // longest-match, so `contig_name` below must try them in this same
    // order and stop at the first success.
    fn m_trinity_comp(b: &[u8]) -> Option<usize> {
        // comp\d+_c\d+_seq\d+
        let p = eat_ci_literal(b, 0, "comp")?;
        let p = eat_digits(b, p)?;
        let p = eat_ci_literal(b, p, "_c")?;
        let p = eat_digits(b, p)?;
        let p = eat_ci_literal(b, p, "_seq")?;
        eat_digits(b, p)
    }

    fn m_trinity_cgi(b: &[u8]) -> Option<usize> {
        // c\d+_g\d+_i\d+
        let p = eat_ci_literal(b, 0, "c")?;
        let p = eat_digits(b, p)?;
        let p = eat_ci_literal(b, p, "_g")?;
        let p = eat_digits(b, p)?;
        let p = eat_ci_literal(b, p, "_i")?;
        eat_digits(b, p)
    }

    fn m_trinity_tr_pipe(b: &[u8]) -> Option<usize> {
        // TR\d+\|c\d+_g\d+_i\d+  (the `\|` is a literal pipe character)
        let p = eat_ci_literal(b, 0, "TR")?;
        let p = eat_digits(b, p)?;
        let p = eat_ci_literal(b, p, "|")?;
        let p = eat_ci_literal(b, p, "c")?;
        let p = eat_digits(b, p)?;
        let p = eat_ci_literal(b, p, "_g")?;
        let p = eat_digits(b, p)?;
        let p = eat_ci_literal(b, p, "_i")?;
        eat_digits(b, p)
    }

    fn m_trinity_dn(b: &[u8]) -> Option<usize> {
        // TRINITY_DN\d+_c\d+_g\d+_i\d+
        let p = eat_ci_literal(b, 0, "TRINITY_DN")?;
        let p = eat_digits(b, p)?;
        let p = eat_ci_literal(b, p, "_c")?;
        let p = eat_digits(b, p)?;
        let p = eat_ci_literal(b, p, "_g")?;
        let p = eat_digits(b, p)?;
        let p = eat_ci_literal(b, p, "_i")?;
        eat_digits(b, p)
    }

    fn m_node(b: &[u8]) -> Option<usize> {
        // node_\d+ (velvet and abyss both use this identical pattern)
        let p = eat_ci_literal(b, 0, "node_")?;
        eat_digits(b, p)
    }

    fn m_idba(b: &[u8]) -> Option<usize> {
        // contig-\d+_\d+
        let p = eat_ci_literal(b, 0, "contig-")?;
        let p = eat_digits(b, p)?;
        let p = eat_ci_literal(b, p, "_")?;
        eat_digits(b, p)
    }

    fn m_spades(b: &[u8]) -> Option<usize> {
        // NODE_\d+_length_\d+_cov_\d+.\d+  (unescaped `.` = any byte)
        let p = eat_ci_literal(b, 0, "NODE_")?;
        let p = eat_digits(b, p)?;
        let p = eat_ci_literal(b, p, "_length_")?;
        let p = eat_digits(b, p)?;
        let p = eat_ci_literal(b, p, "_cov_")?;
        let p = eat_digits(b, p)?;
        let p = eat_any_byte(b, p)?;
        eat_digits(b, p)
    }

    pub fn contig_name(header: &str) -> Option<&str> {
        let b = header.as_bytes();
        let end = m_trinity_comp(b)
            .or_else(|| m_trinity_cgi(b))
            .or_else(|| m_trinity_tr_pipe(b))
            .or_else(|| m_trinity_dn(b))
            .or_else(|| m_node(b))
            .or_else(|| m_idba(b))
            .or_else(|| m_spades(b))?;
        Some(&header[..end])
    }
}

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
    #[error("contig {contig:?} matched {count} UCE loci; expected exactly one")]
    InvalidFilteredMatchCount { contig: String, count: usize },
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
    // ASCII-only guard: every hand-rolled matcher advances by *bytes*,
    // while `regex`'s `.` (used in the header matcher, not this one, but
    // kept consistent here) advances by one Unicode scalar. Restricting
    // the fast path to ASCII input keeps byte offsets == char offsets, so
    // slicing can never panic or land on the wrong boundary.
    if header.is_ascii() && regex.as_str() == fast_extract::DEFAULT_PROBE_REGEX_SRC {
        if let Some(name) = fast_extract::probe_name(header) {
            return Ok(name.to_string());
        }
    }
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
    // See the ASCII-only guard note in `extract_probe_name`; it matters
    // here in particular because the spades pattern's unescaped `.` needs
    // byte-for-byte behavior to match `eat_any_byte`.
    if header.is_ascii() && header_regex.as_str() == fast_extract::default_header_regex_source() {
        if let Some(name) = fast_extract::contig_name(header) {
            return Ok(name.to_string());
        }
    }
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
    pub matches: FastMap<String, FastSet<String>>,
    /// UCE locus -> set of strands it was hit on (usually a single strand)
    pub orientation: FastMap<String, FastSet<String>>,
    /// UCE locus -> set of contigs that hit it
    pub revmatches: FastMap<String, FastSet<String>>,
    /// UCE loci excluded because they're flagged as probe duplicates
    pub probe_dupes: FastSet<String>,
}

/// Mirrors the `main()` per-taxon LASTZ-result loop: classify each match by
/// contig/locus, tracking probe duplicates separately.
pub fn process_taxon_lastz(
    lastz_matches: &[LastzMatch],
    probe_regex: &Regex,
    header_regex: &Regex,
    dupes: &FastSet<String>,
    dupefile_active: bool,
) -> Result<TaxonMatches, MatchError> {
    let mut result = TaxonMatches::default();
    for m in lastz_matches {
        add_taxon_lastz_match(
            &mut result,
            m,
            probe_regex,
            header_regex,
            dupes,
            dupefile_active,
        )?;
    }
    Ok(result)
}

/// Streaming form of [`process_taxon_lastz`], used for large per-taxon
/// alignment files so match rows need not be retained after classification.
pub fn process_taxon_lastz_iter<I>(
    lastz_matches: I,
    probe_regex: &Regex,
    header_regex: &Regex,
    dupes: &FastSet<String>,
    dupefile_active: bool,
) -> Result<TaxonMatches, MatchError>
where
    I: IntoIterator<Item = Result<LastzMatch, phyluce_io::lastz::LastzError>>,
{
    let mut result = TaxonMatches::default();
    for m in lastz_matches {
        let m = m?;
        add_taxon_lastz_match(
            &mut result,
            &m,
            probe_regex,
            header_regex,
            dupes,
            dupefile_active,
        )?;
    }
    Ok(result)
}

fn add_taxon_lastz_match(
    result: &mut TaxonMatches,
    m: &LastzMatch,
    probe_regex: &Regex,
    header_regex: &Regex,
    dupes: &FastSet<String>,
    dupefile_active: bool,
) -> Result<(), MatchError> {
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
    Ok(())
}

/// Contigs that matched more than one distinct UCE locus. Mirrors
/// `check_contigs_for_dupes`.
pub fn contigs_matching_multiple_uces(
    matches: &FastMap<String, FastSet<String>>,
) -> FastSet<String> {
    matches
        .iter()
        .filter(|(_, uces)| uces.len() > 1)
        .map(|(node, _)| node.clone())
        .collect()
}

/// (contigs, UCE loci) where a UCE locus matched more than one contig.
/// Mirrors `check_loci_for_dupes`.
pub fn loci_matching_multiple_contigs(
    revmatches: &FastMap<String, FastSet<String>>,
) -> (FastSet<String>, FastSet<String>) {
    let mut dupe_contigs = FastSet::default();
    let mut dupe_uces = FastSet::default();
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
) -> Result<FastSet<String>, MatchError> {
    let mut matches: FastMap<String, Vec<String>> = FastMap::default();
    for m in lastz_matches {
        let target = extract_probe_name(&m.name1, probe_regex)?;
        let query = extract_probe_name(&m.name2, probe_regex)?;
        matches.entry(target).or_default().push(query);
    }
    let mut dupes = FastSet::default();
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
        matches: &FastMap<String, FastSet<String>>,
        orientation: &FastMap<String, FastSet<String>>,
        critter: &str,
    ) -> Result<(), MatchError> {
        let column = ident(critter);
        let tx = conn.unchecked_transaction()?;
        for (contig, uces) in matches {
            let Some(uce) = uces.iter().next().filter(|_| uces.len() == 1) else {
                return Err(MatchError::InvalidFilteredMatchCount {
                    contig: contig.clone(),
                    count: uces.len(),
                });
            };
            tx.execute(
                &format!("UPDATE matches SET {column} = 1 WHERE uce = ?1"),
                [uce],
            )?;
            let strand = orientation
                .get(uce)
                .and_then(|s| s.iter().next())
                .cloned()
                .unwrap_or_default();
            let orient_value = format!("{contig}({strand})");
            tx.execute(
                &format!("UPDATE match_map SET {column} = ?1 WHERE uce = ?2"),
                rusqlite::params![orient_value, uce],
            )?;
        }
        tx.commit()?;
        Ok(())
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn probe_regex() -> Regex {
        Regex::new(r"^(uce-\d+)(?:_p\d+.*)").unwrap()
    }

    // --- Differential fuzz tests: fast_extract vs the real regex --------
    //
    // The fast path in extract_probe_name/extract_contig_name is only
    // trustworthy if it agrees with `Regex::captures` on every input, not
    // just the happy-path examples above. These generate thousands of
    // inputs -- valid instances of every alternative (varied digit counts,
    // case), truncated/corrupted near-misses, pure garbage, and non-ASCII
    // strings -- and assert byte-for-byte agreement with the regex oracle.
    // A tiny xorshift64 PRNG avoids pulling in a `rand`/`proptest`
    // dependency just for this.

    struct Xorshift64(u64);
    impl Xorshift64 {
        fn new(seed: u64) -> Self {
            Xorshift64(seed | 1)
        }
        fn next_u64(&mut self) -> u64 {
            let mut x = self.0;
            x ^= x << 13;
            x ^= x >> 7;
            x ^= x << 17;
            self.0 = x;
            x
        }
        fn range(&mut self, n: usize) -> usize {
            (self.next_u64() % n as u64) as usize
        }
        fn digits(&mut self, min: usize, max: usize) -> String {
            let len = min + self.range(max - min + 1);
            (0..len)
                .map(|_| (b'0' + self.range(10) as u8) as char)
                .collect()
        }
        fn mutate_case(&mut self, s: &str) -> String {
            s.chars()
                .map(|c| {
                    if self.range(2) == 0 {
                        c.to_ascii_uppercase()
                    } else {
                        c.to_ascii_lowercase()
                    }
                })
                .collect()
        }
        fn garbage(&mut self, len: usize) -> String {
            const CHARS: &[u8] = b"abcXYZ019_-.| ";
            (0..len)
                .map(|_| CHARS[self.range(CHARS.len())] as char)
                .collect()
        }
    }

    fn header_alternative_samples(rng: &mut Xorshift64) -> Vec<String> {
        let mut out = Vec::new();
        // One valid instance of each `[headers]` alternative, several
        // digit-width/case variants each, sometimes with trailing garbage
        // (which the regex's trailing `.*` should still accept).
        for _ in 0..40 {
            let ds: Vec<String> = (0..4).map(|_| rng.digits(1, 4)).collect();
            let variants = [
                format!("comp{}_c{}_seq{}", ds[0], ds[1], ds[2]),
                format!("c{}_g{}_i{}", ds[0], ds[1], ds[2]),
                format!("TR{}|c{}_g{}_i{}", ds[0], ds[1], ds[2], ds[3]),
                format!("TRINITY_DN{}_c{}_g{}_i{}", ds[0], ds[1], ds[2], ds[3]),
                format!("node_{}", ds[0]),
                format!("contig-{}_{}", ds[0], ds[1]),
                format!("NODE_{}_length_{}_cov_{}.{}", ds[0], ds[1], ds[2], ds[3]),
            ];
            for v in variants {
                let v = rng.mutate_case(&v);
                if rng.range(2) == 0 {
                    let garbage_len = rng.range(10);
                    let garbage = rng.garbage(garbage_len);
                    out.push(format!("{v}{garbage}"));
                } else {
                    out.push(v);
                }
            }
        }
        out
    }

    fn probe_samples(rng: &mut Xorshift64) -> Vec<String> {
        let mut out = Vec::new();
        for _ in 0..200 {
            let uce_digits = rng.digits(1, 5);
            let p_digits = rng.digits(1, 3);
            let base = format!("uce-{uce_digits}_p{p_digits}");
            let base = rng.mutate_case(&base);
            let trailing_len = rng.range(15);
            let trailing = rng.garbage(trailing_len);
            out.push(format!("{base}{trailing}"));
        }
        out
    }

    fn corrupted_near_misses(rng: &mut Xorshift64, valid: &[String]) -> Vec<String> {
        // Truncate, drop a char, or insert a char into otherwise-valid
        // samples, to stress the "almost matches" boundary where a
        // hand-rolled scanner is most likely to diverge from the regex.
        let mut out = Vec::new();
        for v in valid {
            if v.is_empty() {
                continue;
            }
            let cut = rng.range(v.len());
            out.push(v[..cut].to_string());
            let del = rng.range(v.len());
            let mut s = v.clone();
            s.remove(del);
            out.push(s);
            let ins = rng.range(v.len() + 1);
            let mut s = v.clone();
            s.insert(ins, (b'a' + rng.range(26) as u8) as char);
            out.push(s);
        }
        out
    }

    fn pure_garbage(rng: &mut Xorshift64) -> Vec<String> {
        let mut out = Vec::with_capacity(500);
        for _ in 0..500 {
            let len = rng.range(30);
            out.push(rng.garbage(len));
        }
        out
    }

    fn non_ascii_samples(rng: &mut Xorshift64) -> Vec<String> {
        // Must never panic (byte-slicing on a non-char-boundary) and must
        // always take the regex path (fast path is ASCII-only).
        vec![
            "uce-42_p1_源".to_string(),
            "NODE_1_length_2_cov_3.é4".to_string(),
            "comp1_c2_seq3🧬".to_string(),
            format!("uce-{}_p1", rng.digits(1, 3)) + "\u{0301}",
        ]
    }

    const FUZZ_SEEDS: &[u64] = &[
        0xC0FFEE, 0xDEADBEEF, 1, 2, 42, 0x1234_5678_9ABC_DEF0, 0xFFFF_FFFF, 0x5EED_5EED,
    ];

    #[test]
    fn fast_extract_probe_name_matches_regex_oracle() {
        let re = probe_regex();
        assert_eq!(re.as_str(), fast_extract::DEFAULT_PROBE_REGEX_SRC);

        for &seed in FUZZ_SEEDS {
            let mut rng = Xorshift64::new(seed);
            let mut samples = probe_samples(&mut rng);
            samples.extend(corrupted_near_misses(&mut rng, &samples.clone()));
            samples.extend(pure_garbage(&mut rng));
            samples.extend(non_ascii_samples(&mut rng));

            for header in &samples {
                let expected = re
                    .captures(header)
                    .and_then(|c| c.get(1))
                    .map(|m| m.as_str().to_string());
                let actual = extract_probe_name(header, &re).ok();
                assert_eq!(expected, actual, "diverged on header={header:?} (seed={seed:#x})");
            }
        }
    }

    #[test]
    fn fast_extract_contig_name_matches_regex_oracle() {
        let src = fast_extract::default_header_regex_source();
        assert!(!src.is_empty(), "packaged default config failed to load");
        let re = Regex::new(src).unwrap();

        for &seed in FUZZ_SEEDS {
            let mut rng = Xorshift64::new(seed);
            let mut samples = header_alternative_samples(&mut rng);
            samples.extend(corrupted_near_misses(&mut rng, &samples.clone()));
            samples.extend(pure_garbage(&mut rng));
            samples.extend(non_ascii_samples(&mut rng));

            for header in &samples {
                let expected = re
                    .captures(header)
                    .and_then(|c| c.get(1))
                    .map(|m| m.as_str().to_string());
                let actual = extract_contig_name(header, &re).ok();
                assert_eq!(expected, actual, "diverged on header={header:?} (seed={seed:#x})");
            }
        }
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

    // --- Ad hoc benchmarks (not part of the normal test run) ------------
    //
    // Exploring where `match-contigs-to-probes` actually spends time before
    // committing to a hashing-library swap. Run with:
    //   cargo +stable test --release -p phyluce-assembly --lib -- --ignored --nocapture bench_
    //
    // Synthetic workload: ~300k LASTZ match rows against ~5k contigs and
    // ~1.5k UCE loci, sized to resemble a mid-size non-model-organism
    // capture run (many contigs, each hitting a handful of loci).

    fn synthetic_lastz_rows(n: usize) -> Vec<phyluce_io::lastz::LastzMatch> {
        (0..n)
            .map(|i| {
                let contig = i % 5_000;
                let uce = i % 1_500;
                phyluce_io::lastz::LastzMatch {
                    score: "4000".to_string(),
                    name1: format!("NODE_{contig}_length_500_cov_10.3"),
                    strand1: "+".to_string(),
                    zstart1: 0,
                    end1: 100,
                    length1: 500,
                    name2: format!("uce-{uce}_p1 |source:test"),
                    strand2: if i % 2 == 0 { "+" } else { "-" }.to_string(),
                    zstart2: 0,
                    end2: 100,
                    length2: 100,
                    diff: String::new(),
                    cigar: "100M".to_string(),
                    identity: "95.0".to_string(),
                    percent_identity: 95.0,
                    continuity: "100.0".to_string(),
                    percent_continuity: 100.0,
                    coverage: None,
                    percent_coverage: None,
                }
            })
            .collect()
    }

    #[test]
    #[ignore]
    fn bench_process_taxon_lastz_end_to_end() {
        let rows = synthetic_lastz_rows(300_000);
        let probe_re = probe_regex();
        // A byte-identical-to-default header regex, so the fast_extract
        // path actually engages -- unlike the plain `node_\d+` used
        // elsewhere in these tests, which is a real, common config but not
        // the packaged default, so it deliberately takes the regex-only
        // path.
        let header_re = Regex::new(fast_extract::default_header_regex_source()).unwrap();
        let dupes = FastSet::default();

        let start = std::time::Instant::now();
        let result = process_taxon_lastz(&rows, &probe_re, &header_re, &dupes, false).unwrap();
        let elapsed = start.elapsed();

        eprintln!(
            "[bench] process_taxon_lastz (ahash + fast_extract): {} rows -> {} contigs matched in {:?} ({:.0} rows/sec)",
            rows.len(),
            result.matches.len(),
            elapsed,
            rows.len() as f64 / elapsed.as_secs_f64()
        );
    }

    #[test]
    #[ignore]
    fn bench_hashmap_std_vs_fxhash_isolated() {
        // Isolates just the HashMap<String, HashSet<String>>
        // entry().or_default().insert() pattern `add_taxon_lastz_match`
        // does per row, with the regex/formatting cost already stripped
        // out, to see how much of the total is attributable to hashing.
        let n = 300_000usize;
        let contig_keys: Vec<String> = (0..n).map(|i| format!("NODE_{}", i % 5_000)).collect();
        let uce_vals: Vec<String> = (0..n).map(|i| format!("uce-{}", i % 1_500)).collect();

        let start = std::time::Instant::now();
        let mut std_map: HashMap<String, HashSet<String>> = HashMap::new();
        for i in 0..n {
            std_map
                .entry(contig_keys[i].clone())
                .or_default()
                .insert(uce_vals[i].clone());
        }
        let std_elapsed = start.elapsed();

        let start = std::time::Instant::now();
        let mut fx_map: rustc_hash::FxHashMap<String, rustc_hash::FxHashSet<String>> =
            rustc_hash::FxHashMap::default();
        for i in 0..n {
            fx_map
                .entry(contig_keys[i].clone())
                .or_default()
                .insert(uce_vals[i].clone());
        }
        let fx_elapsed = start.elapsed();

        let start = std::time::Instant::now();
        let mut ah_map: FastMap<String, FastSet<String>> = FastMap::default();
        for i in 0..n {
            ah_map
                .entry(contig_keys[i].clone())
                .or_default()
                .insert(uce_vals[i].clone());
        }
        let ah_elapsed = start.elapsed();

        assert_eq!(std_map.len(), fx_map.len());
        assert_eq!(std_map.len(), ah_map.len());
        eprintln!(
            "[bench] {n} inserts: std HashMap {:?} vs FxHashMap {:?} ({:.2}x) vs AHashMap {:?} ({:.2}x)",
            std_elapsed,
            fx_elapsed,
            std_elapsed.as_secs_f64() / fx_elapsed.as_secs_f64(),
            ah_elapsed,
            std_elapsed.as_secs_f64() / ah_elapsed.as_secs_f64(),
        );
    }

    #[test]
    #[ignore]
    fn bench_regex_captures_vs_captures_read() {
        // Isolates whether `Regex::captures` (allocates a fresh `Captures`
        // per call) vs `Regex::captures_read` (reuses one `CaptureLocations`
        // buffer across calls) matters here, as a lower-risk alternative to
        // hand-rolling a byte scanner: both `probe_regex` (`--regex`) and
        // `header_regex` (built from phyluce.conf's `[headers]` section,
        // a multi-way case-insensitive alternation across assemblers) are
        // user-configurable, so a hand-rolled fast path would only cover
        // the default patterns and still need a regex fallback for custom
        // ones. `captures_read` keeps the general regex engine but drops
        // the per-call allocation.
        let rows = synthetic_lastz_rows(300_000);
        let probe_re = probe_regex();
        let header_re = contig_header_regex(r"node_\d+").unwrap();

        let start = std::time::Instant::now();
        let mut total_len = 0usize;
        for row in &rows {
            let c1 = header_re.captures(&row.name1).unwrap();
            let contig = c1.get(1).unwrap().as_str();
            let c2 = probe_re.captures(&row.name2).unwrap();
            let uce = c2.get(1).unwrap().as_str();
            total_len += contig.len() + uce.len();
        }
        let captures_elapsed = start.elapsed();

        let mut header_locs = header_re.capture_locations();
        let mut probe_locs = probe_re.capture_locations();
        let start = std::time::Instant::now();
        let mut total_len2 = 0usize;
        for row in &rows {
            header_re.captures_read(&mut header_locs, &row.name1).unwrap();
            let (s1, e1) = header_locs.get(1).unwrap();
            let contig = &row.name1[s1..e1];
            probe_re.captures_read(&mut probe_locs, &row.name2).unwrap();
            let (s2, e2) = probe_locs.get(1).unwrap();
            let uce = &row.name2[s2..e2];
            total_len2 += contig.len() + uce.len();
        }
        let captures_read_elapsed = start.elapsed();

        assert_eq!(total_len, total_len2);
        eprintln!(
            "[bench] {} rows: captures() {:?} vs captures_read() {:?} ({:.2}x)",
            rows.len(),
            captures_elapsed,
            captures_read_elapsed,
            captures_elapsed.as_secs_f64() / captures_read_elapsed.as_secs_f64()
        );
    }

    #[test]
    #[ignore]
    fn bench_regex_vs_hand_rolled_default_pattern_ceiling() {
        // Quantifies the *ceiling*: a hand-rolled scanner that only
        // understands the single default `^(uce-\d+)(?:_p\d+.*)` pattern
        // (no support for a user's `--regex` override), vs the general
        // regex engine. Not meant to be shipped as-is -- see
        // bench_regex_captures_vs_captures_read's doc comment for why a
        // hand-rolled path would need a regex fallback anyway.
        fn hand_rolled_uce_prefix(header: &str) -> Option<&str> {
            let rest = header.strip_prefix("uce-")?;
            let digit_len = rest.bytes().take_while(u8::is_ascii_digit).count();
            if digit_len == 0 {
                return None;
            }
            let after = &rest[digit_len..];
            if after.starts_with("_p") && after[2..].bytes().next()?.is_ascii_digit() {
                Some(&header[..4 + digit_len])
            } else {
                None
            }
        }

        let rows = synthetic_lastz_rows(300_000);
        let probe_re = probe_regex();

        let start = std::time::Instant::now();
        let mut total_len = 0usize;
        for row in &rows {
            let c = probe_re.captures(&row.name2).unwrap();
            total_len += c.get(1).unwrap().as_str().len();
        }
        let regex_elapsed = start.elapsed();

        let start = std::time::Instant::now();
        let mut total_len2 = 0usize;
        for row in &rows {
            total_len2 += hand_rolled_uce_prefix(&row.name2).unwrap().len();
        }
        let hand_rolled_elapsed = start.elapsed();

        assert_eq!(total_len, total_len2);
        eprintln!(
            "[bench] {} rows, probe extraction only: regex {:?} vs hand-rolled {:?} ({:.2}x)",
            rows.len(),
            regex_elapsed,
            hand_rolled_elapsed,
            regex_elapsed.as_secs_f64() / hand_rolled_elapsed.as_secs_f64()
        );
    }

    #[test]
    #[ignore]
    fn bench_regex_extraction_isolated() {
        // Isolates just the two `Regex::captures` calls `add_taxon_lastz_match`
        // does per row (extract_contig_name / extract_probe_name), to see
        // how much of the ~168ms full-pipeline time (see
        // bench_process_taxon_lastz_end_to_end) is regex vs hashing.
        let rows = synthetic_lastz_rows(300_000);
        let probe_re = probe_regex();
        let header_re = contig_header_regex(r"node_\d+").unwrap();

        let start = std::time::Instant::now();
        let mut total_len = 0usize;
        for row in &rows {
            let contig = extract_contig_name(&row.name1, &header_re).unwrap();
            let uce = extract_probe_name(&row.name2, &probe_re).unwrap();
            total_len += contig.len() + uce.len();
        }
        let elapsed = start.elapsed();

        eprintln!(
            "[bench] regex extraction only: {} rows in {:?} ({:.0} rows/sec, checksum={total_len})",
            rows.len(),
            elapsed,
            rows.len() as f64 / elapsed.as_secs_f64()
        );
    }

    #[test]
    fn detects_contigs_matching_multiple_uces() {
        let mut matches = FastMap::default();
        matches.insert(
            "NODE_1".to_string(),
            FastSet::from_iter(["uce-1".to_string(), "uce-2".to_string()]),
        );
        matches.insert("NODE_2".to_string(), FastSet::from_iter(["uce-3".to_string()]));
        let dupes = contigs_matching_multiple_uces(&matches);
        assert_eq!(dupes, FastSet::from_iter(["NODE_1".to_string()]));
    }

    #[test]
    fn detects_loci_matching_multiple_contigs() {
        let mut revmatches = FastMap::default();
        revmatches.insert(
            "uce-1".to_string(),
            FastSet::from_iter(["NODE_1".to_string(), "NODE_2".to_string()]),
        );
        revmatches.insert("uce-2".to_string(), FastSet::from_iter(["NODE_3".to_string()]));
        let (dupe_contigs, dupe_uces) = loci_matching_multiple_contigs(&revmatches);
        assert_eq!(
            dupe_contigs,
            FastSet::from_iter(["NODE_1".to_string(), "NODE_2".to_string()])
        );
        assert_eq!(dupe_uces, FastSet::from_iter(["uce-1".to_string()]));
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

        let mut matches = FastMap::default();
        matches.insert("NODE_1".to_string(), FastSet::from_iter(["uce-1".to_string()]));
        let mut orientation = FastMap::default();
        orientation.insert("uce-1".to_string(), FastSet::from_iter(["+".to_string()]));
        db::store_lastz_results(&conn, &matches, &orientation, "taxon_a").unwrap();

        let invalid: FastMap<String, FastSet<String>> = FastMap::from_iter([(
            "NODE_2".to_string(),
            FastSet::from_iter(["uce-1".to_string(), "uce-2".to_string()]),
        )]);
        assert!(matches!(
            db::store_lastz_results(&conn, &invalid, &orientation, "taxon_b"),
            Err(MatchError::InvalidFilteredMatchCount { count: 2, .. })
        ));

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
