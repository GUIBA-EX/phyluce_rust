//! CLI wiring for `phyluce probe get-multi-fasta-table` /
//! `phyluce probe query-multi-fasta-table`, mirroring
//! `phyluce_probe_get_multi_fasta_table` / `phyluce_probe_query_multi_fasta_table`.

use std::collections::HashMap;
use std::path::{Path, PathBuf};

use phyluce_io::read_fasta;
use phyluce_io::sql::ident;
use rusqlite::Connection;

/// Mirrors `create_match_database` + the insert loop.
pub fn run_get(fastas_dir: &Path, output: &Path, base_taxon: &str) -> anyhow::Result<()> {
    let mut fasta_files: Vec<PathBuf> = std::fs::read_dir(fastas_dir)?
        .filter_map(|e| e.ok())
        .map(|e| e.path())
        .filter(|p| p.extension().and_then(|e| e.to_str()) == Some("fasta"))
        .collect();
    fasta_files.sort();

    let mut organisms = Vec::new();
    let mut conserved: HashMap<String, Vec<String>> = HashMap::new();
    println!("Reading Fasta files...");
    for fasta in &fasta_files {
        let taxon_name = fasta
            .file_stem()
            .and_then(|s| s.to_str())
            .unwrap_or("")
            .to_string();
        organisms.push(taxon_name.clone());
        let records = read_fasta(fasta)?;
        for record in &records {
            // mirrors `seq.id.split("|")[3].split(":")[1]`
            let locus = record
                .id
                .split('|')
                .nth(3)
                .and_then(|s| s.split_once(':'))
                .map(|(_, v)| v.to_string())
                .ok_or_else(|| anyhow::anyhow!("record '{}': missing locus metadata", record.id))?;
            conserved.entry(locus).or_default().push(taxon_name.clone());
        }
    }

    if output.exists() {
        std::fs::remove_file(output)?;
    }
    let conn = Connection::open(output)?;
    conn.execute_batch("PRAGMA foreign_keys = ON;")?;
    let table = ident(base_taxon);
    let columns = organisms
        .iter()
        .map(|o| format!("{} integer DEFAULT 0", ident(o)))
        .collect::<Vec<_>>()
        .join(", ");
    conn.execute(
        &format!("CREATE TABLE {table} (locus text primary key, {columns})"),
        [],
    )?;

    for (locus, taxa) in &conserved {
        let names = taxa.iter().map(|t| ident(t)).collect::<Vec<_>>().join(", ");
        let ones = vec!["1"; taxa.len()].join(", ");
        conn.execute(
            &format!("INSERT INTO {table} (locus, {names}) values (?1, {ones})"),
            [locus],
        )?;
    }
    Ok(())
}

/// Mirrors the query script's default (non `--specific-counts`) summary
/// path, plus the `--specific-counts` filtering path.
pub fn run_query(
    db: &Path,
    base_taxon: &str,
    specific_counts: Option<usize>,
    output: Option<&Path>,
) -> anyhow::Result<()> {
    if let Some(_sc) = specific_counts {
        anyhow::ensure!(output.is_some(), "--specific-counts requires --output");
    }
    let conn = Connection::open(db)?;
    let table = ident(base_taxon);
    let mut stmt = conn.prepare(&format!("PRAGMA table_info({table})"))?;
    let columns: Vec<String> = stmt
        .query_map([], |r| r.get::<_, String>(1))?
        .collect::<Result<_, _>>()?;
    let taxa = &columns[1..];

    let mut stmt = conn.prepare(&format!("SELECT * FROM {table}"))?;
    let column_count = columns.len();

    if let Some(threshold) = specific_counts {
        let threshold = threshold as i64;
        let output = output.unwrap();
        let missing_path = format!("{}.missing.matrix", output.display());
        let mut out1 = std::fs::File::create(output)?;
        let mut out2 = std::fs::File::create(&missing_path)?;
        use std::io::Write as _;
        writeln!(out1, "# Hits against {threshold} taxa")?;
        writeln!(out1, "# {}", taxa.join(","))?;
        writeln!(out1, "[hits]")?;
        writeln!(out2, "# {}", taxa.join(","))?;
        writeln!(out2, "[misses]")?;

        let mut counter: HashMap<String, usize> = HashMap::new();
        let mut total_loci = 0usize;
        let mut rows = stmt.query([])?;
        while let Some(row) = rows.next()? {
            let locus: String = row.get(0)?;
            let mut values = Vec::new();
            let mut sum = 0i64;
            for i in 1..column_count {
                let v: i64 = row.get(i)?;
                values.push(v);
                sum += v;
            }
            if sum >= threshold {
                for (i, &v) in values.iter().enumerate() {
                    if v == 1 {
                        *counter.entry(taxa[i].clone()).or_insert(0) += 1;
                    }
                }
                total_loci += 1;
                writeln!(out1, "{locus}")?;
            } else {
                let joined = values
                    .iter()
                    .map(|v| v.to_string())
                    .collect::<Vec<_>>()
                    .join(",");
                writeln!(out2, "{locus}\t{sum}\t{joined}")?;
            }
        }
        println!("{counter:?}");
        println!("Total loci = {total_loci}");
    } else {
        let mut counts = vec![0usize; taxa.len() + 1];
        let mut rows = stmt.query([])?;
        while let Some(row) = rows.next()? {
            let mut sum = 0i64;
            for i in 1..column_count {
                sum += row.get::<_, i64>(i)?;
            }
            counts[sum as usize] += 1;
        }
        for i in 0..=taxa.len() {
            let total: usize = counts[i..].iter().sum();
            println!("Loci shared by {i} taxa:\t{total}");
        }
    }
    Ok(())
}
