//! CLI wiring for `phyluce align get-taxon-locus-counts-in-alignments`,
//! mirroring `phyluce_align_get_taxon_locus_counts_in_alignments`.

use std::collections::HashMap;
use std::path::Path;

use crate::informative_sites_cmd::{find_alignment_files, load_alignment};

pub fn run(alignments_dir: &Path, input_format: &str, output: &Path) -> anyhow::Result<()> {
    let files = find_alignment_files(alignments_dir, input_format)?;

    let mut order: Vec<String> = Vec::new();
    let mut counts: HashMap<String, usize> = HashMap::new();
    for file in &files {
        let alignment = load_alignment(file, input_format)?;
        for row in &alignment.rows {
            let distinct: std::collections::HashSet<u8> = row.seq.iter().copied().collect();
            if distinct.len() > 1 {
                if !counts.contains_key(&row.id) {
                    order.push(row.id.clone());
                }
                *counts.entry(row.id.clone()).or_insert(0) += 1;
            }
        }
    }

    crate::cli_info!("Writing taxon count data to {}", output.display());
    let mut out = std::fs::File::create(output)?;
    use std::io::Write as _;
    writeln!(out, "taxon,count")?;
    for taxon in &order {
        writeln!(out, "{taxon},{}", counts[taxon])?;
    }
    Ok(())
}
