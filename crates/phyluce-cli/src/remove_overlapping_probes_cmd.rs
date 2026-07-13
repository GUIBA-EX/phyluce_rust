//! CLI wiring for `phyluce probe remove-overlapping-probes-given-config`,
//! mirroring `phyluce_probe_remove_overlapping_probes_given_config`.

use std::collections::HashSet;
use std::path::Path;

use phyluce_io::{read_fasta, write_fasta_record};

pub fn run(probes: &Path, config: &Path, output: &Path) -> anyhow::Result<()> {
    let conf_text = std::fs::read_to_string(config)?;
    let sections = crate::conf::parse_ini(&conf_text);
    let excludes: HashSet<String> = sections
        .get("exclude")
        .map(|entries| entries.iter().map(|(k, _)| k.clone()).collect())
        .unwrap_or_default();
    crate::cli_info!("There are {} loci to exclude", excludes.len());

    let records = read_fasta(probes)?;
    let mut kept_loci = HashSet::new();
    let mut dropped_loci = HashSet::new();
    let mut out = std::fs::File::create(output)?;
    for record in &records {
        let locus = record
            .id
            .split('_')
            .next()
            .unwrap_or(&record.id)
            .to_string();
        if !excludes.contains(&locus) {
            write_fasta_record(&mut out, &record.description, &record.sequence)?;
            kept_loci.insert(locus);
        } else {
            dropped_loci.insert(locus);
        }
    }
    crate::cli_info!("Kept {} loci", kept_loci.len());
    crate::cli_info!("Dropped {} loci", dropped_loci.len());
    Ok(())
}
