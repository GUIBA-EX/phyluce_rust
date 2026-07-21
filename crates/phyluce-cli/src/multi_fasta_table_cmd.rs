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
    crate::cli_info!("Reading Fasta files...");
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

    // One transaction for the whole insert loop instead of autocommit's
    // implicit per-statement transaction/fsync: 5000 single-row INSERTs
    // without this take ~3.3s, ~4.7ms with it (~700x) -- see
    // `tests::bench_sqlite_insert_autocommit_vs_one_transaction`.
    let tx = conn.unchecked_transaction()?;
    for (locus, taxa) in &conserved {
        let names = taxa.iter().map(|t| ident(t)).collect::<Vec<_>>().join(", ");
        let ones = vec!["1"; taxa.len()].join(", ");
        tx.execute(
            &format!("INSERT INTO {table} (locus, {names}) values (?1, {ones})"),
            [locus],
        )?;
    }
    tx.commit()?;
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
    anyhow::ensure!(
        columns.len() >= 2 && columns.first().is_some_and(|name| name == "locus"),
        "database table {base_taxon:?} must contain a 'locus' column and at least one taxon column"
    );
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
        crate::cli_info!("{counter:?}");
        crate::cli_info!("Total loci = {total_loci}");
    } else {
        let mut counts = vec![0usize; taxa.len() + 1];
        let mut rows = stmt.query([])?;
        while let Some(row) = rows.next()? {
            let mut sum = 0i64;
            for i in 1..column_count {
                sum += row.get::<_, i64>(i)?;
            }
            anyhow::ensure!(
                (0..=taxa.len() as i64).contains(&sum),
                "row in {base_taxon:?} has taxon-column sum {sum}, expected a value in 0..={} (are all taxon columns 0 or 1?)",
                taxa.len()
            );
            counts[sum as usize] += 1;
        }
        for i in 0..=taxa.len() {
            let total: usize = counts[i..].iter().sum();
            crate::cli_info!("Loci shared by {i} taxa:\t{total}");
        }
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    // Ad hoc benchmark: `run_get`'s insert loop calls `conn.execute(INSERT
    // ...)` once per locus with no explicit transaction, so each INSERT is
    // its own SQLite autocommit transaction -- one fsync per row on disk,
    // regardless of how fast the in-process algorithm is. Compares against
    // wrapping the same inserts in one transaction, matching the pattern
    // `phyluce-assembly::db::store_lastz_results` already uses elsewhere.
    // Uses a real file-backed connection (not `:memory:`), since
    // `:memory:` never touches disk and would hide exactly the cost this
    // is measuring. Run with:
    //   cargo +stable test --release -p phyluce-cli --bin phyluce -- --ignored --nocapture bench_sqlite_insert
    use rusqlite::Connection;

    fn temp_db_path(name: &str) -> std::path::PathBuf {
        let dir = std::env::temp_dir().join("phyluce-cli-sqlite-bench");
        std::fs::create_dir_all(&dir).unwrap();
        dir.join(format!("{name}-{}.sqlite", std::process::id()))
    }

    #[test]
    #[ignore]
    fn bench_sqlite_insert_autocommit_vs_one_transaction() {
        let n = 5_000;

        let path = temp_db_path("autocommit");
        let _ = std::fs::remove_file(&path);
        let conn = Connection::open(&path).unwrap();
        conn.execute("CREATE TABLE t (id INTEGER PRIMARY KEY, v TEXT)", [])
            .unwrap();
        let start = std::time::Instant::now();
        for i in 0..n {
            conn.execute(
                "INSERT INTO t (id, v) VALUES (?1, ?2)",
                rusqlite::params![i, "x"],
            )
            .unwrap();
        }
        let autocommit_elapsed = start.elapsed();
        drop(conn);
        std::fs::remove_file(&path).ok();

        let path = temp_db_path("one-tx");
        let _ = std::fs::remove_file(&path);
        let conn = Connection::open(&path).unwrap();
        conn.execute("CREATE TABLE t (id INTEGER PRIMARY KEY, v TEXT)", [])
            .unwrap();
        let start = std::time::Instant::now();
        let tx = conn.unchecked_transaction().unwrap();
        for i in 0..n {
            tx.execute(
                "INSERT INTO t (id, v) VALUES (?1, ?2)",
                rusqlite::params![i, "x"],
            )
            .unwrap();
        }
        tx.commit().unwrap();
        let tx_elapsed = start.elapsed();
        drop(conn);
        std::fs::remove_file(&path).ok();

        eprintln!(
            "[bench] {n} inserts: autocommit (no tx) {:?} vs one transaction {:?} ({:.1}x)",
            autocommit_elapsed,
            tx_elapsed,
            autocommit_elapsed.as_secs_f64() / tx_elapsed.as_secs_f64()
        );
    }
}
