//! CLI wiring for `phyluce utilities filter-bed-by-fasta`, mirroring
//! `phyluce_utilities_filter_bed_by_fasta`.

use std::collections::HashSet;
use std::io::Write as _;
use std::path::{Path, PathBuf};

use phyluce_io::read_fasta;

pub fn run(bed: &Path, fasta: &Path, output: Option<PathBuf>) -> anyhow::Result<()> {
    let records = read_fasta(fasta)?;
    let loci: HashSet<String> = records
        .iter()
        .filter_map(|r| r.description.split('|').nth(1))
        .map(|s| s.trim().to_string())
        .collect();

    let bed_text = std::fs::read_to_string(bed)?;
    let mut out: Box<dyn std::io::Write> = match &output {
        Some(p) => Box::new(std::fs::File::create(p)?),
        None => Box::new(std::io::stdout()),
    };
    for line in bed_text.lines() {
        if line.starts_with("track") {
            writeln!(out, "{line}")?;
            continue;
        }
        let fields: Vec<&str> = line.split('\t').collect();
        if let Some(name_field) = fields.get(3) {
            let locus = name_field.split('_').next().unwrap_or("");
            if loci.contains(locus) {
                writeln!(out, "{line}")?;
            }
        }
    }
    Ok(())
}
