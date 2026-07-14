//! CLI wiring for `phyluce align get-align-summary-data`, mirroring
//! `phyluce_align_get_align_summary_data`'s `--output-stats` CSV (the log
//! output isn't reproduced byte-for-byte; it isn't covered by any golden
//! fixture).

use std::path::{Path, PathBuf};

use phyluce_align::summary::compute_align_summary;

use crate::informative_sites_cmd::{find_alignment_files, load_alignment};

pub fn run(
    alignments_dir: &Path,
    input_format: &str,
    output_stats: Option<PathBuf>,
    show_taxon_counts: bool,
    cores: usize,
) -> anyhow::Result<()> {
    anyhow::ensure!(cores > 0, "--cores must be greater than zero");
    let files = find_alignment_files(alignments_dir, input_format)?;
    anyhow::ensure!(!files.is_empty(), "no {input_format} alignments found");
    let mut rows = Vec::new();
    for file in &files {
        let name = file
            .file_name()
            .and_then(|n| n.to_str())
            .unwrap_or("")
            .to_string();
        let alignment = load_alignment(file, input_format)?;
        let s = compute_align_summary(&alignment);
        let missing = alignment
            .rows
            .iter()
            .map(|row| row.seq.iter().filter(|&&c| c == b'?').count())
            .sum::<usize>() as f64
            / (alignment.ntax() * alignment.nchar()) as f64;
        rows.push((name, alignment.ntax(), missing, s));
    }

    print_aggregate_summary(&rows, show_taxon_counts)?;

    if let Some(out_path) = output_stats {
        let mut out = String::from(
            "aln,length,sites,differences,characters,gc content,gaps,a count, c count, g count, t count\n",
        );
        for (name, _, _, s) in &rows {
            out.push_str(&format!(
                "{},{},{},{},{},{},{},{},{},{},{}\n",
                name,
                s.length,
                s.sum_informative_sites,
                s.sum_differences,
                s.sum_counted_sites,
                format_gc(s.gc_content_percent()),
                s.char_count(b'-'),
                s.char_count(b'A'),
                s.char_count(b'C'),
                s.char_count(b'G'),
                s.char_count(b'T'),
            ));
        }
        std::fs::write(out_path, out)?;
    }
    Ok(())
}

fn mean(values: &[f64]) -> f64 {
    values.iter().sum::<f64>() / values.len() as f64
}

fn ci95(values: &[f64]) -> f64 {
    if values.len() <= 1 {
        return f64::NAN;
    }
    let average = mean(values);
    let variance = values
        .iter()
        .map(|value| (value - average).powi(2))
        .sum::<f64>()
        / (values.len() - 1) as f64;
    1.96 * variance.sqrt() / (values.len() as f64).sqrt()
}

fn print_aggregate_summary(
    rows: &[(String, usize, f64, phyluce_align::summary::AlignSummary)],
    show_taxon_counts: bool,
) -> anyhow::Result<()> {
    let lengths: Vec<f64> = rows.iter().map(|(_, _, _, s)| s.length as f64).collect();
    let sites: Vec<f64> = rows
        .iter()
        .map(|(_, _, _, s)| s.sum_informative_sites as f64)
        .collect();
    let taxa: Vec<f64> = rows.iter().map(|(_, n, _, _)| *n as f64).collect();
    let missing: Vec<f64> = rows.iter().map(|(_, _, m, _)| *m).collect();

    crate::cli_info!("[Alignments] loci:\t{}", rows.len());
    crate::cli_info!(
        "[Alignments] length:\t{}",
        lengths.iter().sum::<f64>() as usize
    );
    crate::cli_info!("[Alignments] mean:\t{:.2}", mean(&lengths));
    crate::cli_info!("[Alignments] 95% CI:\t{:.2}", ci95(&lengths));
    crate::cli_info!(
        "[Alignments] min:\t{}",
        lengths.iter().copied().fold(f64::INFINITY, f64::min) as usize
    );
    crate::cli_info!(
        "[Alignments] max:\t{}",
        lengths.iter().copied().fold(f64::NEG_INFINITY, f64::max) as usize
    );
    crate::cli_info!("[Sites] total:\t{}", sites.iter().sum::<f64>() as usize);
    crate::cli_info!("[Sites] mean:\t{:.2}", mean(&sites));
    crate::cli_info!("[Sites] 95% CI:\t{:.2}", ci95(&sites));
    crate::cli_info!(
        "[Sites] min:\t{}",
        sites.iter().copied().fold(f64::INFINITY, f64::min) as usize
    );
    crate::cli_info!(
        "[Sites] max:\t{}",
        sites.iter().copied().fold(f64::NEG_INFINITY, f64::max) as usize
    );
    crate::cli_info!("[Taxa] mean:\t{:.2}", mean(&taxa));
    crate::cli_info!("[Taxa] 95% CI:\t{:.2}", ci95(&taxa));
    crate::cli_info!(
        "[Taxa] min:\t{}",
        taxa.iter().copied().fold(f64::INFINITY, f64::min) as usize
    );
    crate::cli_info!(
        "[Taxa] max:\t{}",
        taxa.iter().copied().fold(f64::NEG_INFINITY, f64::max) as usize
    );
    crate::cli_info!("[Missing] mean:\t{:.2}", mean(&missing) * 100.0);
    crate::cli_info!("[Missing] 95% CI:\t{:.2}", ci95(&missing) * 100.0);
    crate::cli_info!(
        "[Missing] min:\t{:.2}",
        missing.iter().copied().fold(f64::INFINITY, f64::min) * 100.0
    );
    crate::cli_info!(
        "[Missing] max:\t{:.2}",
        missing.iter().copied().fold(f64::NEG_INFINITY, f64::max) * 100.0
    );

    let mut characters = std::collections::BTreeMap::new();
    for (_, _, _, summary) in rows {
        for (&character, &count) in &summary.characters {
            *characters.entry(character).or_insert(0usize) += count;
        }
    }
    let all_characters: usize = characters.values().sum();
    let nucleotides: usize = [b'A', b'C', b'G', b'T']
        .iter()
        .map(|base| characters.get(base).copied().unwrap_or(0))
        .sum();
    crate::cli_info!("[All characters]\t{all_characters}");
    crate::cli_info!("[Nucleotides]\t\t{nucleotides}");
    for (character, count) in characters {
        crate::cli_info!(
            "[Characters] '{}' is present {count} times",
            character as char
        );
    }

    let max_taxa = taxa.iter().copied().fold(0.0, f64::max) as usize;
    for step in 10..20 {
        let percent = step as f64 * 0.05;
        let minimum = ((percent - 0.01) * max_taxa as f64).ceil() as usize;
        let count = rows.iter().filter(|(_, n, _, _)| *n >= minimum).count();
        crate::cli_info!(
            "[Matrix {}%]\t\t{count} alignments",
            (percent * 100.0) as usize
        );
    }

    if show_taxon_counts {
        let mut counts = std::collections::BTreeMap::new();
        for (_, ntax, _, _) in rows {
            *counts.entry(*ntax).or_insert(0usize) += 1;
        }
        for (ntax, count) in counts {
            crate::cli_info!("[Taxa] {count} alignments contain {ntax} taxa");
        }
    }
    Ok(())
}

/// Mirrors Python's default float `str()` formatting used by `"{}".format`
/// on a `round(x, 2)` result: whole numbers print without a trailing `.0`
/// stripped (e.g. `50.0`, not `50`), matching `round()`'s float return type.
fn format_gc(x: f64) -> String {
    if x == x.trunc() {
        format!("{x:.1}")
    } else {
        format!("{x}")
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn summary_helpers_match_numpy_sample_statistics() {
        let values = [10.0, 20.0, 30.0];
        assert_eq!(mean(&values), 20.0);
        assert!((ci95(&values) - 11.316_065).abs() < 0.000_001);
        assert!(ci95(&[10.0]).is_nan());
    }
}
