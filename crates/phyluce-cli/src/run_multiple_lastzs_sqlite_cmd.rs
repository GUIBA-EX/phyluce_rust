//! CLI wiring for `phyluce probe run-multiple-lastzs-sqlite`, mirroring
//! `phyluce_probe_run_multiple_lastzs_sqlite`.
//!
//! Untested: `lastz` isn't installed in this environment. Also simplified
//! relative to the Python original: `phyluce.many_lastz.multi_lastz_runner`
//! chunks/parallelizes the target genome (per-sequence for chromosomes,
//! size-bucketed temp FASTAs for "huge" scaffold sets) across
//! `multiprocessing.Pool(cores)` workers and concatenates their outputs.
//! This port runs a single `lastz` invocation per genome directly against
//! the whole `.2bit` file instead -- `--cores` is accepted but unused. The
//! resulting `.lastz` file and SQLite tables have the same shape either
//! way; only the parallel chunking strategy (irrelevant to correctness,
//! only to wall-clock time on huge genomes) is not reproduced.

use std::path::{Path, PathBuf};

use phyluce_io::sql::ident;
use rusqlite::Connection;

pub struct RunMultipleLastzsArgs {
    pub chromolist: Vec<String>,
    pub scaffoldlist: Vec<String>,
    pub append: bool,
    pub no_dir: bool,
    pub genome_base_path: String,
    pub coverage: f64,
    pub identity: f64,
}

fn genome_path(base_path: &str, no_dir: bool, name: &str) -> PathBuf {
    if no_dir {
        Path::new(base_path).join(format!("{name}.2bit"))
    } else {
        Path::new(base_path).join(name).join(format!("{name}.2bit"))
    }
}

fn create_species_lastz_table(conn: &Connection, g: &str) -> anyhow::Result<()> {
    let table = ident(g);
    let index = ident(&format!("{g}_name2_idx"));
    conn.execute_batch(&format!(
        "CREATE TABLE {table} (
            id INTEGER PRIMARY KEY AUTOINCREMENT,
            score INTEGER NOT NULL,
            name1 TEXT NOT NULL,
            strand1 TEXT NOT NULL,
            zstart1 INTEGER NOT NULL,
            end1 INTEGER NOT NULL,
            length1 INTEGER NOT NULL,
            name2 TEXT NOT NULL,
            strand2 TEXT NOT NULL,
            zstart2 INTEGER NOT NULL,
            end2 INTEGER NOT NULL,
            length2 INTEGER NOT NULL,
            diff TEXT NOT NULL,
            cigar TEXT NOT NULL,
            identity TEXT NOT NULL,
            percent_identity FLOAT NOT NULL,
            continuity TEXT NOT NULL,
            percent_continuity FLOAT NOT NULL,
            coverage TEXT NOT NULL,
            percent_coverage FLOAT NOT NULL);
         CREATE INDEX {index} on {table}(name2);"
    ))?;
    Ok(())
}

/// Mirrors `clean_lastz_data`: strip literal `%` characters (lastz writes
/// percentages like `100.0%` in several columns; the DB stores them as
/// bare numbers).
fn clean_lastz_data(text: &str) -> String {
    text.replace('%', "")
}

fn insert_species_rows(conn: &Connection, g: &str, cleaned: &str) -> anyhow::Result<()> {
    let table = ident(g);
    let mut stmt = conn.prepare(&format!(
        "INSERT INTO {table} (score, name1, strand1, zstart1, end1,
            length1, name2, strand2, zstart2, end2, length2, diff, cigar,
            identity, percent_identity, continuity, percent_continuity,
            coverage, percent_coverage) VALUES (?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?)"
    ))?;
    for line in cleaned.lines() {
        let fields: Vec<&str> = line.trim().split('\t').collect();
        anyhow::ensure!(fields.len() == 19, "malformed lastz row: {line:?}");
        stmt.execute(rusqlite::params_from_iter(fields.iter()))?;
    }
    Ok(())
}

fn store_species_rows(conn: &Connection, genome: &str, cleaned: &str) -> anyhow::Result<()> {
    let tx = conn.unchecked_transaction()?;
    create_species_lastz_table(&tx, genome)?;
    insert_species_rows(&tx, genome, cleaned)?;
    tx.execute("INSERT INTO species (name) VALUES (?1)", [genome])?;
    tx.commit()?;
    Ok(())
}

fn align_and_store(
    conn: &Connection,
    lastz_bin: &str,
    genomes: &[String],
    probefile: &Path,
    output_dir: &Path,
    args: &RunMultipleLastzsArgs,
) -> anyhow::Result<()> {
    for g in genomes {
        let target = genome_path(&args.genome_base_path, args.no_dir, g);
        let probe_name = probefile
            .file_name()
            .and_then(|s| s.to_str())
            .unwrap_or("probes");
        let output =
            crate::output_path::output_file(output_dir, &format!("{probe_name}_v_{g}.lastz"))?;
        crate::lastz_align::run_many_lastz(
            lastz_bin,
            &target.to_string_lossy(),
            &probefile.to_string_lossy(),
            args.coverage,
            args.identity,
            &output.to_string_lossy(),
        )?;
        let raw = std::fs::read_to_string(&output)?;
        let cleaned = clean_lastz_data(&raw);
        let clean_path = format!("{}.clean", output.display());
        std::fs::write(&clean_path, &cleaned)?;
        std::fs::remove_file(&output)?;

        store_species_rows(conn, g, &cleaned)?;
    }
    Ok(())
}

pub fn run(
    db: &Path,
    output_dir: &Path,
    probefile: &Path,
    args: &RunMultipleLastzsArgs,
) -> anyhow::Result<()> {
    std::fs::create_dir_all(output_dir)?;
    for g in args.chromolist.iter().chain(args.scaffoldlist.iter()) {
        let path = genome_path(&args.genome_base_path, args.no_dir, g);
        anyhow::ensure!(path.is_file(), "{} is not a file", path.display());
    }

    let cfg = phyluce_config::PhyluceConfig::load()?;
    let lastz_bin = cfg.get_user_path("binaries", "lastz")?;

    let conn = Connection::open(db)?;
    if !args.append {
        conn.execute(
            "CREATE TABLE species (name TEXT PRIMARY KEY, description TEXT NULL, version TEXT NULL)",
            [],
        )?;
    }

    align_and_store(
        &conn,
        &lastz_bin,
        &args.scaffoldlist,
        probefile,
        output_dir,
        args,
    )?;
    align_and_store(
        &conn,
        &lastz_bin,
        &args.chromolist,
        probefile,
        output_dir,
        args,
    )?;
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn clean_lastz_data_strips_percent_signs() {
        assert_eq!(clean_lastz_data("100.0%\t83.2%\tfoo"), "100.0\t83.2\tfoo");
    }

    #[test]
    fn genome_path_respects_no_dir() {
        assert_eq!(
            genome_path("/genomes", false, "gallus"),
            PathBuf::from("/genomes/gallus/gallus.2bit")
        );
        assert_eq!(
            genome_path("/genomes", true, "gallus"),
            PathBuf::from("/genomes/gallus.2bit")
        );
    }

    #[test]
    fn species_import_rolls_back_on_a_malformed_row() {
        let conn = Connection::open_in_memory().unwrap();
        conn.execute("CREATE TABLE species (name TEXT PRIMARY KEY)", [])
            .unwrap();
        let valid = vec!["1"; 19].join("\t");
        let input = format!("{valid}\nmalformed\n");
        assert!(store_species_rows(&conn, "taxon_a", &input).is_err());

        let table_count: i64 = conn
            .query_row(
                "SELECT count(*) FROM sqlite_master WHERE type='table' AND name='taxon_a'",
                [],
                |row| row.get(0),
            )
            .unwrap();
        let species_count: i64 = conn
            .query_row("SELECT count(*) FROM species", [], |row| row.get(0))
            .unwrap();
        assert_eq!(table_count, 0);
        assert_eq!(species_count, 0);
    }
}
