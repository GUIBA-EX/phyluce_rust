//! CLI wiring for `phyluce align get-incomplete-matrix-estimates`,
//! mirroring `phyluce_align_get_incomplete_matrix_estimates`.

use std::collections::BTreeMap;
use std::path::Path;

use phyluce_io::sql::ident;
use rusqlite::Connection;

fn validate_range(min: f64, max: f64, step: f64) -> anyhow::Result<()> {
    anyhow::ensure!(
        min.is_finite() && (0.0..1.0).contains(&min),
        "The min value must be 0 <= value < 1"
    );
    anyhow::ensure!(
        max.is_finite() && (0.0..=1.0).contains(&max),
        "The max value must be 0 <= value <= 1"
    );
    anyhow::ensure!(
        step.is_finite() && step > 0.0 && step <= 1.0,
        "The step value must be 0 < value <= 1"
    );
    anyhow::ensure!(min < max, "The min value must be less than max");
    Ok(())
}

pub fn run(
    db: &Path,
    min: f64,
    max: f64,
    step: f64,
    exclude: &[String],
    include: &[String],
) -> anyhow::Result<()> {
    validate_range(min, max, step)?;

    let conn = Connection::open(db)?;
    let mut stmt = conn.prepare("PRAGMA table_info(matches)")?;
    let all_columns: Vec<String> = stmt
        .query_map([], |row| row.get::<_, String>(1))?
        .collect::<Result<_, _>>()?;
    anyhow::ensure!(
        all_columns.len() >= 2 && all_columns.first().is_some_and(|name| name == "uce"),
        "database table 'matches' must contain an 'uce' column and at least one taxon column"
    );
    // first column is the `uce` id column, matching Python's `all_columns[1:]`
    let candidate_taxa = &all_columns[1..];

    let taxa: Vec<String> = if !exclude.is_empty() {
        let excludes: std::collections::HashSet<&String> = exclude.iter().collect();
        candidate_taxa
            .iter()
            .filter(|t| !excludes.contains(t))
            .cloned()
            .collect()
    } else if !include.is_empty() {
        let includes: std::collections::HashSet<&String> = include.iter().collect();
        candidate_taxa
            .iter()
            .filter(|t| includes.contains(t))
            .cloned()
            .collect()
    } else {
        candidate_taxa.to_vec()
    };
    crate::cli_info!("There are {} taxa.", taxa.len());

    let cols = taxa.iter().map(|t| ident(t)).collect::<Vec<_>>().join(",");
    let query = format!("SELECT uce,{cols} FROM matches");
    let mut stmt = conn.prepare(&query)?;
    let n = taxa.len();
    let mut locus_counts: BTreeMap<String, i64> = BTreeMap::new();
    let mut rows = stmt.query([])?;
    while let Some(row) = rows.next()? {
        let locus: String = row.get(0)?;
        let mut sum = 0i64;
        for i in 0..n {
            let v: String = row.get(i + 1)?;
            if v == "1" {
                sum += 1;
            }
        }
        locus_counts.insert(locus, sum);
    }

    let num_taxa = taxa.len() as f64;
    let mut fracs = Vec::new();
    let mut f = min;
    while f < max {
        fracs.push(f);
        f += step;
    }
    let cuts: Vec<f64> = fracs.iter().map(|f| (num_taxa * f).round()).collect();

    let mut frac_dict: BTreeMap<String, usize> = BTreeMap::new();
    for count in locus_counts.values() {
        let v = *count as f64;
        let idx = cuts.iter().position(|&c| v < c);
        let key = match idx {
            None => format!(
                "{},{}",
                format_frac(*fracs.last().unwrap()),
                format_cut(*cuts.last().unwrap())
            ),
            Some(0) => {
                // Python: `positions[0] - 1` with positions[0] == 0 wraps to
                // index -1 (last element) via negative indexing; reproduced
                // verbatim even though it looks like an off-by-one bug.
                format!(
                    "{},{}",
                    format_frac(*fracs.last().unwrap()),
                    format_cut(*cuts.last().unwrap())
                )
            }
            Some(p) => format!("{},{}", format_frac(fracs[p - 1]), format_cut(cuts[p - 1])),
        };
        *frac_dict.entry(key).or_insert(0) += 1;
    }

    crate::cli_info!("Freq(taxa present),Cut point,Loci");
    for (k, v) in &frac_dict {
        crate::cli_info!("{k},{v}");
    }
    Ok(())
}

fn format_frac(f: f64) -> String {
    // Python's repr of numpy floats from `arange`/`around` prints the
    // shortest round-trippable decimal; Rust's default float Display is a
    // reasonable approximation for the typical 1-2 decimal steps phyluce uses.
    format!("{f}")
}

fn format_cut(c: f64) -> String {
    format!("{c}")
}

#[cfg(test)]
mod tests {
    use super::{run, validate_range};

    #[test]
    fn rejects_non_progressing_and_empty_ranges() {
        assert!(validate_range(0.5, 0.9, 0.0).is_err());
        assert!(validate_range(0.5, 0.9, -0.1).is_err());
        assert!(validate_range(0.9, 0.5, 0.1).is_err());
        assert!(validate_range(0.5, 0.5, 0.1).is_err());
        assert!(validate_range(f64::NAN, 0.9, 0.1).is_err());
        assert!(validate_range(0.5, 0.9, 0.1).is_ok());
    }

    #[test]
    fn rejects_database_without_matches_schema() {
        let path = std::env::temp_dir().join(format!(
            "phyluce-empty-matrix-{}.sqlite",
            std::process::id()
        ));
        let _ = std::fs::remove_file(&path);
        assert!(run(&path, 0.0, 1.0, 0.1, &[], &[]).is_err());
        std::fs::remove_file(path).unwrap();
    }
}
