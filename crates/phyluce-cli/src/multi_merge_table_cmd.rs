//! CLI wiring for `phyluce probe get-multi-merge-table` /
//! `phyluce probe query-multi-merge-table`, mirroring
//! `phyluce_probe_get_multi_merge_table` / `phyluce_probe_query_multi_merge_table`.

use std::collections::HashMap;
use std::path::Path;

use phyluce_io::sql::ident;
use rusqlite::Connection;

struct Interval {
    start: i64,
    end: i64,
    taxa: std::collections::HashSet<String>,
}

/// Mirrors the `bx.intervals.intersection.IntervalTree`-based accumulation:
/// a new interval either joins every *existing* interval on the same
/// chromosome that it overlaps (adding its taxon to each), or starts a new
/// interval if there's no overlap.
fn insert_interval(
    tree: &mut HashMap<String, Vec<Interval>>,
    chromo: &str,
    start: i64,
    end: i64,
    taxon: &str,
) {
    let intervals = tree.entry(chromo.to_string()).or_default();
    let mut any_overlap = false;
    for iv in intervals.iter_mut() {
        if start < iv.end && iv.start < end {
            iv.taxa.insert(taxon.to_string());
            any_overlap = true;
        }
    }
    if !any_overlap {
        intervals.push(Interval {
            start,
            end,
            taxa: std::collections::HashSet::from([taxon.to_string()]),
        });
    }
}

pub fn run_get(conf: &Path, output: &Path, base_taxon: &str) -> anyhow::Result<()> {
    let conf_text = std::fs::read_to_string(conf)?;
    let sections = crate::conf::parse_ini(&conf_text);
    let beds = sections
        .get("beds")
        .ok_or_else(|| anyhow::anyhow!("no [beds] section in --conf"))?;

    let mut organisms = Vec::new();
    let mut conserved: HashMap<String, Vec<Interval>> = HashMap::new();
    crate::cli_info!("Reading the BED file for:");
    for (taxon_name, bedfile) in beds {
        organisms.push(taxon_name.clone());
        let text = std::fs::read_to_string(bedfile)?;
        for line in text.lines() {
            let fields: Vec<&str> = line.trim().split('\t').collect();
            if fields.len() < 3 {
                continue;
            }
            let chromo = fields[0];
            let start: i64 = fields[1].parse()?;
            let stop: i64 = fields[2].parse()?;
            insert_interval(&mut conserved, chromo, start, stop, taxon_name);
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
        &format!(
            "CREATE TABLE {table} (uce integer primary key autoincrement, chromo text, start integer, stop integer, {columns})"
        ),
        [],
    )?;

    let mut chromos: Vec<&String> = conserved.keys().collect();
    chromos.sort();
    // One transaction for the whole insert loop -- see the identical fix
    // (and benchmark) in `multi_fasta_table_cmd::run_get`: per-statement
    // autocommit means one fsync per row, ~700x slower than a single
    // transaction at a few thousand rows.
    let tx = conn.unchecked_transaction()?;
    for chromo in chromos {
        let mut intervals: Vec<&Interval> = conserved[chromo].iter().collect();
        intervals.sort_by_key(|iv| iv.start);
        for iv in intervals {
            let mut names: Vec<&String> = iv.taxa.iter().collect();
            names.sort();
            let names_joined = names
                .iter()
                .map(|s| ident(s))
                .collect::<Vec<_>>()
                .join(", ");
            let ones = vec!["1"; names.len()].join(", ");
            tx.execute(
                &format!(
                    "INSERT INTO {table}(chromo, start, stop, {names_joined}) values (?1, ?2, ?3, {ones})"
                ),
                rusqlite::params![chromo, iv.start, iv.end],
            )?;
        }
    }
    tx.commit()?;
    Ok(())
}

pub fn run_query(
    db: &Path,
    base_taxon: &str,
    specific_counts: Option<usize>,
    output: Option<&Path>,
) -> anyhow::Result<()> {
    if specific_counts.is_some() {
        anyhow::ensure!(output.is_some(), "--specific-counts requires --output");
    }
    let conn = Connection::open(db)?;
    let table = ident(base_taxon);
    let mut stmt = conn.prepare(&format!("PRAGMA table_info({table})"))?;
    let columns: Vec<String> = stmt
        .query_map([], |r| r.get::<_, String>(1))?
        .collect::<Result<_, _>>()?;
    // columns: uce, chromo, start, stop, <taxa...>
    anyhow::ensure!(
        columns.len() >= 5
            && columns[..4] == ["uce", "chromo", "start", "stop"],
        "database table {base_taxon:?} must contain uce/chromo/start/stop and at least one taxon column"
    );
    let taxa = &columns[4..];
    let column_count = columns.len();

    let mut stmt = conn.prepare(&format!("SELECT * FROM {table}"))?;

    if let Some(threshold) = specific_counts {
        let threshold = threshold as i64;
        let output = output.unwrap();
        let missing_path = format!("{}.missing.matrix", output.display());
        let mut out1 = std::fs::File::create(output)?;
        let mut out2 = std::fs::File::create(&missing_path)?;
        use std::io::Write as _;
        writeln!(out2, "{}", taxa.join(","))?;

        let mut counter: HashMap<String, usize> = HashMap::new();
        let mut rows = stmt.query([])?;
        while let Some(row) = rows.next()? {
            let chromo: String = row.get(1)?;
            let start: i64 = row.get(2)?;
            let stop: i64 = row.get(3)?;
            let mut values = Vec::new();
            let mut sum = 0i64;
            for i in 4..column_count {
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
                writeln!(out1, "{chromo}\t{start}\t{stop}")?;
            } else {
                let joined = values
                    .iter()
                    .map(|v| v.to_string())
                    .collect::<Vec<_>>()
                    .join(",");
                writeln!(out2, "{joined}")?;
            }
        }
        crate::cli_info!("{counter:?}");
    } else {
        let mut counts = vec![0usize; taxa.len() + 1];
        let mut rows = stmt.query([])?;
        while let Some(row) = rows.next()? {
            let mut sum = 0i64;
            for i in 4..column_count {
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
            crate::cli_info!("Loci shared by {base_taxon} + {i} taxa:\t{total}");
        }
    }
    Ok(())
}
