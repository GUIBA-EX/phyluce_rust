//! CLI wiring for `phyluce probe get-subsets-of-tiled-probes`, mirroring
//! `phyluce_probe_get_subsets_of_tiled_probes`.

use std::collections::HashMap;
use std::path::Path;

use anyhow::Context;
use phyluce_io::{read_fasta, write_fasta_record};
use regex::Regex;

pub fn run(probes: &Path, taxa: &[String], output: &Path, regex_str: &str) -> anyhow::Result<()> {
    let regex = Regex::new(regex_str)?;
    let records =
        read_fasta(probes).with_context(|| format!("reading probes {}", probes.display()))?;

    let mut taxa_counts: HashMap<String, usize> = HashMap::new();
    let mut kept_loci: std::collections::HashSet<String> = std::collections::HashSet::new();
    let mut probes_kept = 0usize;

    let mut out = std::fs::File::create(output)
        .with_context(|| format!("creating output file {}", output.display()))?;
    for record in &records {
        // mirrors `seq.description.split("|")[1].split(",")[4].split(":")[1]`
        let field =
            record.description.split('|').nth(1).ok_or_else(|| {
                anyhow::anyhow!("record '{}': missing '|' metadata field", record.id)
            })?;
        let taxon = field
            .split(',')
            .nth(4)
            .and_then(|s| s.split_once(':'))
            .map(|(_, v)| v.to_string())
            .ok_or_else(|| anyhow::anyhow!("record '{}': missing taxon metadata", record.id))?;
        *taxa_counts.entry(taxon.clone()).or_insert(0) += 1;

        if taxa.contains(&taxon) {
            write_fasta_record(&mut out, &record.description, &record.sequence)?;
            probes_kept += 1;
            let locus = regex
                .captures(&record.id)
                .and_then(|c| c.get(1))
                .map(|m| m.as_str().to_string())
                .ok_or_else(|| anyhow::anyhow!("no regex match for probe id {:?}", record.id))?;
            kept_loci.insert(locus);
        }
    }

    crate::cli_info!("All probes = {}", taxa_counts.values().sum::<usize>());
    crate::cli_info!("--- Probes by taxon ---");
    for (taxon, count) in &taxa_counts {
        crate::cli_info!("{taxon}\t{count}");
    }
    crate::cli_info!("--- Post  filtering ---");
    crate::cli_info!("Conserved locus count = {}", kept_loci.len());
    crate::cli_info!("Probe Count = {probes_kept}");
    Ok(())
}
