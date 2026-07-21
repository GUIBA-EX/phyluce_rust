//! CLI wiring for `phyluce align get-smilogram-from-alignments`, mirroring
//! `phyluce_align_get_smilogram_from_alignments`.
//!
//! The Python original breaks a genuine allele-count tie via
//! `random.choice`; this port always picks the first-encountered
//! candidate instead -- a deliberate, documented divergence, not a bug.
//! It also silently overwrites `--output-database` if it already exists,
//! instead of Python's interactive "Overwrite [Y/n]?" prompt (which can't
//! usefully run non-interactively anyway).

use std::collections::HashMap;
use std::io::Write as _;
use std::path::Path;

use anyhow::Context;
use rusqlite::Connection;

use crate::informative_sites_cmd::{find_alignment_files, load_alignment};

const IUPAC_CANDIDATES: &[char] = &[
    'a', 'c', 't', 'g', 'r', 'y', 's', 'w', 'k', 'm', 'b', 'd', 'h', 'v',
];

fn replace_gaps_at_ends(seq: &[u8]) -> Vec<u8> {
    let mut out = seq.to_vec();
    let mut start = 0;
    while start < out.len() && out[start] == b'-' {
        out[start] = b'?';
        start += 1;
    }
    let mut end = out.len();
    while end > start && out[end - 1] == b'-' {
        end -= 1;
        out[end] = b'?';
    }
    out
}

#[derive(Default)]
struct TaxonBuckets {
    insertion: Vec<i64>,
    deletion: Vec<i64>,
    substitution: Vec<i64>,
    majallele: Vec<i64>,
    missing: Vec<i64>,
}

fn create_tables(conn: &Connection) -> anyhow::Result<()> {
    conn.execute_batch(
        "CREATE TABLE loci (locus text PRIMARY KEY, length int);
         CREATE TABLE by_taxon (
            idx INTEGER PRIMARY KEY AUTOINCREMENT,
            taxon text, locus text, position int, position_from_center real, type text,
            FOREIGN KEY (locus) REFERENCES loci(locus));
         CREATE TABLE by_locus (
            idx INTEGER PRIMARY KEY AUTOINCREMENT,
            locus text, majallele real, substitutions real, deletions real, insertions real,
            missing real, bases real, position int, position_from_center real, type text,
            FOREIGN KEY (locus) REFERENCES loci(locus));
         CREATE TABLE by_locus_missing (
            idx INTEGER PRIMARY KEY AUTOINCREMENT,
            locus text, present real, absent real, position int, position_from_center real, type text,
            FOREIGN KEY (locus) REFERENCES loci(locus));",
    )?;
    Ok(())
}

#[allow(clippy::type_complexity)]
fn process_alignment(
    seqs: &[(String, Vec<u8>)],
) -> (HashMap<String, TaxonBuckets>, usize, HashMap<usize, usize>) {
    let aligned: Vec<Vec<u8>> = seqs.iter().map(|(_, s)| replace_gaps_at_ends(s)).collect();
    let length = aligned.first().map(|s| s.len()).unwrap_or(0);

    let mut results: HashMap<String, TaxonBuckets> = seqs
        .iter()
        .map(|(id, _)| (id.clone(), TaxonBuckets::default()))
        .collect();
    let mut base_count: HashMap<usize, usize> = HashMap::new();

    for idx in 0..length {
        let col: Vec<u8> = aligned
            .iter()
            .map(|s| s[idx].to_ascii_lowercase())
            .collect();
        let bases: Vec<u8> = col
            .iter()
            .copied()
            .filter(|&b| b != b'n' && b != b'?')
            .collect();
        base_count.insert(idx, bases.len());

        if bases.is_empty() {
            for (row, (id, _)) in seqs.iter().enumerate() {
                let base = col[row];
                if base == b'n' || base == b'?' {
                    results.get_mut(id).unwrap().missing.push(idx as i64);
                }
            }
            continue;
        }
        let major_base = {
            let unique: std::collections::HashSet<u8> = bases.iter().copied().collect();
            if unique.len() == 1 {
                bases[0]
            } else {
                let mut order: Vec<u8> = Vec::new();
                let mut counts: HashMap<u8, usize> = HashMap::new();
                for &b in &bases {
                    if !counts.contains_key(&b) {
                        order.push(b);
                    }
                    *counts.entry(b).or_insert(0) += 1;
                }
                let max_count = *counts.values().max().unwrap();
                let tied: Vec<u8> = order
                    .iter()
                    .copied()
                    .filter(|b| counts[b] == max_count)
                    .collect();
                if tied.len() == 1 {
                    tied[0]
                } else {
                    let candidates: Vec<u8> = order
                        .iter()
                        .copied()
                        .filter(|b| {
                            counts[b] == max_count && IUPAC_CANDIDATES.contains(&(*b as char))
                        })
                        .collect();
                    candidates.first().copied().unwrap_or(tied[0])
                }
            }
        };

        for (row, (id, _)) in seqs.iter().enumerate() {
            let base = col[row];
            let bucket = results.get_mut(id).unwrap();
            if base == b'n' || base == b'?' {
                bucket.missing.push(idx as i64);
            } else if base == major_base {
                bucket.majallele.push(idx as i64);
            } else if major_base == b'-' && base != b'-' {
                bucket.insertion.push(idx as i64);
            } else if base == b'-' && major_base != b'-' {
                bucket.deletion.push(idx as i64);
            } else if base != b'-' && major_base != b'-' {
                bucket.substitution.push(idx as i64);
            }
        }
    }
    (results, length, base_count)
}

pub fn run(
    alignments_dir: &Path,
    output_file: &Path,
    output_missing: &Path,
    output_database: &Path,
    input_format: &str,
) -> anyhow::Result<()> {
    if output_database.exists() {
        std::fs::remove_file(output_database)
            .with_context(|| format!("removing existing database {}", output_database.display()))?;
    }
    let conn = Connection::open(output_database)
        .with_context(|| format!("opening database {}", output_database.display()))?;
    create_tables(&conn)?;

    let files = find_alignment_files(alignments_dir, input_format)
        .with_context(|| format!("reading alignments directory {}", alignments_dir.display()))?;
    // One transaction for every insert below instead of autocommit's
    // implicit per-statement transaction/fsync -- this loop can do a
    // per-(taxon, position) INSERT, so a real dataset (thousands of loci x
    // tens of taxa x several variant positions each) can easily reach
    // hundreds of thousands of individual statements. At the ~700x-slower
    // rate measured for autocommit in `multi_fasta_table_cmd::tests::
    // bench_sqlite_insert_autocommit_vs_one_transaction`, that's the
    // difference between sub-second and (extrapolated) many minutes.
    let tx = conn.unchecked_transaction()?;
    for file in &files {
        let name = file.file_name().and_then(|n| n.to_str()).unwrap_or("");
        let locus = name.split('.').next().unwrap_or(name).to_string();
        let alignment = load_alignment(file, input_format)
            .with_context(|| format!("loading alignment {}", file.display()))?;
        let seqs: Vec<(String, Vec<u8>)> = alignment
            .rows
            .iter()
            .map(|r| (r.id.clone(), r.seq.clone()))
            .collect();
        let (results, length, base_count) = process_alignment(&seqs);
        let center = length as f64 / 2.0;

        tx.execute(
            "INSERT INTO loci VALUES (?1, ?2)",
            rusqlite::params![locus, length as i64],
        )?;

        for (taxon, buckets) in &results {
            for (typ, positions) in [
                (Some("insertion"), &buckets.insertion),
                (Some("deletion"), &buckets.deletion),
                (Some("substitution"), &buckets.substitution),
                (Some("majallele"), &buckets.majallele),
                (None, &buckets.missing),
            ] {
                for &pos in positions {
                    tx.execute(
                        "INSERT INTO by_taxon (taxon, locus, position, position_from_center, type) VALUES (?1,?2,?3,?4,?5)",
                        rusqlite::params![taxon, locus, pos, pos as f64 - center, typ],
                    )?;
                }
            }
        }

        let mut maj_cnt: HashMap<i64, i64> = HashMap::new();
        let mut subs_cnt: HashMap<i64, i64> = HashMap::new();
        let mut dels_cnt: HashMap<i64, i64> = HashMap::new();
        let mut ins_cnt: HashMap<i64, i64> = HashMap::new();
        let mut n_cnt: HashMap<i64, i64> = HashMap::new();
        for buckets in results.values() {
            for &p in &buckets.majallele {
                *maj_cnt.entry(p).or_insert(0) += 1;
            }
            for &p in &buckets.substitution {
                *subs_cnt.entry(p).or_insert(0) += 1;
            }
            for &p in &buckets.deletion {
                *dels_cnt.entry(p).or_insert(0) += 1;
            }
            for &p in &buckets.insertion {
                *ins_cnt.entry(p).or_insert(0) += 1;
            }
            for &p in &buckets.missing {
                *n_cnt.entry(p).or_insert(0) += 1;
            }
        }

        let mut positions: Vec<usize> = base_count.keys().copied().collect();
        positions.sort();
        for pos in positions {
            let p = pos as i64;
            let bases = base_count[&pos] as f64;
            tx.execute(
                "INSERT INTO by_locus (locus, majallele, substitutions, deletions, insertions, missing, bases, position, position_from_center, type)
                 VALUES (?1,?2,?3,?4,?5,?6,?7,?8,?9,?10)",
                rusqlite::params![
                    locus,
                    *maj_cnt.get(&p).unwrap_or(&0) as f64,
                    *subs_cnt.get(&p).unwrap_or(&0) as f64,
                    *dels_cnt.get(&p).unwrap_or(&0) as f64,
                    *ins_cnt.get(&p).unwrap_or(&0) as f64,
                    *n_cnt.get(&p).unwrap_or(&0) as f64,
                    bases,
                    p,
                    p as f64 - center,
                    "substitutions",
                ],
            )?;
            tx.execute(
                "INSERT INTO by_locus_missing (locus, present, absent, position, position_from_center, type) VALUES (?1,?2,?3,?4,?5,?6)",
                rusqlite::params![locus, bases, results.len() as f64 - bases, p, p as f64 - center, "missing"],
            )?;
        }
    }
    tx.commit()?;

    let mut outf = std::fs::File::create(output_file)
        .with_context(|| format!("creating output file {}", output_file.display()))?;
    writeln!(outf, "substitutions,bp,freq,distance_from_center")?;
    conn.execute_batch(
        "CREATE TEMP TABLE ssb AS
         SELECT sum(substitutions) AS ss, sum(bases) AS total_bases, sum(substitutions)/sum(bases) AS freq, position_from_center
         FROM by_locus GROUP BY position_from_center",
    )?;
    {
        let mut stmt = conn
            .prepare("SELECT ss, total_bases, freq, position_from_center FROM ssb WHERE ss != 0")?;
        let mut rows = stmt.query([])?;
        while let Some(row) = rows.next()? {
            let a: f64 = row.get(0)?;
            let b: f64 = row.get(1)?;
            let c: f64 = row.get(2)?;
            let d: f64 = row.get(3)?;
            writeln!(outf, "{a},{b},{c},{d}")?;
        }
    }

    let mut outm = std::fs::File::create(output_missing)
        .with_context(|| format!("creating output file {}", output_missing.display()))?;
    writeln!(outm, "substitutions,bp,freq,distance_from_center")?;
    conn.execute_batch(
        "CREATE TEMP TABLE ssc AS
         SELECT sum(present) as pres, sum(absent) AS total_absent, sum(absent)/(sum(absent) + sum(present)) AS freq, position_from_center
         FROM by_locus_missing GROUP BY position_from_center",
    )?;
    {
        let mut stmt = conn.prepare(
            "SELECT pres, total_absent, freq, position_from_center FROM ssc WHERE pres != 0",
        )?;
        let mut rows = stmt.query([])?;
        while let Some(row) = rows.next()? {
            let a: f64 = row.get(0)?;
            let b: f64 = row.get(1)?;
            let c: f64 = row.get(2)?;
            let d: f64 = row.get(3)?;
            writeln!(outm, "{a},{b},{c},{d}")?;
        }
    }
    Ok(())
}
