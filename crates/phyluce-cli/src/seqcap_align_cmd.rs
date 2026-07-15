//! CLI wiring for `phyluce align seqcap-align`, mirroring
//! `phyluce_align_seqcap_align` (MAFFT path; `--aligner muscle` isn't
//! implemented yet).

use std::collections::HashMap;
use std::path::Path;

use anyhow::Context;
use phyluce_align::mafft::{mafft_binary_path, run_mafft};
use phyluce_align::nexus::format_nexus;
use phyluce_align::trim::{trim_alignment_running, validate_trim_parameters};
use phyluce_align::Alignment;
use phyluce_config::PhyluceConfig;
use phyluce_io::read_fasta;

#[allow(clippy::too_many_arguments)]
pub fn run(
    input: &Path,
    output: &Path,
    taxa: usize,
    incomplete_matrix: bool,
    no_trim: bool,
    ambiguous: bool,
    window: usize,
    proportion: f64,
    threshold: f64,
    max_divergence: f64,
    min_length: usize,
    cores: usize,
) -> anyhow::Result<()> {
    anyhow::ensure!(cores > 0, "--cores must be greater than zero");
    validate_trim_parameters(window, proportion, threshold, max_divergence)?;
    crate::output_path::prepare_output_dir(output)?;

    if ambiguous {
        crate::cli_info!("NOT removing sequences with ambiguous bases...");
    } else {
        crate::cli_info!("Removing ALL sequences with ambiguous bases...");
    }

    let records = read_fasta(input)?;
    let mut loci: HashMap<String, Vec<(String, String)>> = HashMap::new();
    for record in records {
        // mirrors `record.description.split("|")[1].rstrip("_phased")`
        let locus = record
            .description
            .split('|')
            .nth(1)
            .map(|s| s.trim_end_matches("_phased").to_string())
            .with_context(|| format!("record '{}' has no '|locus' field", record.id))?;
        let has_n = record.sequence.contains('N') || record.sequence.contains('n');
        if ambiguous || !has_n {
            loci.entry(locus)
                .or_default()
                .push((record.id, record.sequence));
        } else {
            crate::cli_warn!("Skipping {locus} because it contains ambiguous bases");
        }
    }

    let min_taxa = if incomplete_matrix { 3 } else { taxa };
    let mut locus_names: Vec<String> = loci.keys().cloned().collect();
    locus_names.sort();
    for locus in &locus_names {
        if loci[locus].len() < min_taxa {
            loci.remove(locus);
            if incomplete_matrix {
                crate::cli_warn!("DROPPED locus {locus}. Too few taxa (N < 3).");
            } else {
                crate::cli_warn!(
                    "DROPPED locus {locus}. Alignment does not contain all {taxa} taxa."
                );
            }
        }
    }

    let cfg = PhyluceConfig::load()?;
    let mafft_bin = mafft_binary_path(&cfg)?;

    let mut work: Vec<(String, Vec<(String, String)>)> = loci.into_iter().collect();
    work.sort_by(|left, right| left.0.cmp(&right.0));
    let results = crate::parallel::try_map_ordered(work, cores, |(locus, sequences)| {
        let raw = run_mafft(&mafft_bin, &sequences)?;
        let aligned = if no_trim {
            Some(raw)
        } else {
            trim_alignment_running(
                &raw,
                window,
                proportion,
                threshold,
                max_divergence,
                min_length,
            )
        };
        match aligned {
            Some(aln) => {
                write_output(output, &locus, &aln)?;
                Ok(false)
            }
            None => Ok(true),
        }
    })?;
    for &dropped in &results {
        print!("{}", if dropped { 'X' } else { '.' });
    }
    crate::cli_info!();
    let dropped = results.iter().filter(|&&dropped| dropped).count();
    if dropped > 0 {
        crate::cli_info!("Dropped {dropped} alignment(s)");
    }
    Ok(())
}

fn write_output(output_dir: &Path, locus: &str, alignment: &Alignment) -> anyhow::Result<()> {
    let out_path = crate::output_path::output_file(output_dir, &format!("{locus}.nexus"))?;
    std::fs::write(out_path, format_nexus(alignment))?;
    Ok(())
}
