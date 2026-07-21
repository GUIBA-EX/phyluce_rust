//! CLI wiring for `phyluce probe reconstruct-uce-from-probe`, mirroring
//! `phyluce_probe_reconstruct_uce_from_probe`.
//!
//! Multi-probe loci use MAFFT by default and retain the Python original's
//! Biopython `dumb_consensus()` semantics. MUSCLE 3 `-clwstrict` remains
//! available as an explicit legacy compatibility path.
//!
//! `--per-locus-dir` is an addition beyond the Python original: alongside
//! the usual single combined `--output` FASTA, it also writes one
//! unwrapped, single-record FASTA file per locus (`<locus>.fasta`) into
//! the given directory -- the reference-directory layout GeneMiner2-UCE's
//! `-r` flag expects (one file per UCE, p1/p2/... overlap already merged
//! into a single sequence by the same consensus step used for `--output`).

use std::collections::HashMap;
use std::io::Write as _;
use std::path::Path;
use std::process::{Command, Stdio};

use anyhow::Context;
use phyluce_io::FastaRecord;

/// Mirrors `Bio.Align.AlignInfo.SummaryInfo.dumb_consensus` with its
/// default `threshold=0.7`, `ambiguous="X"`.
fn dumb_consensus(seqs: &[Vec<u8>]) -> Vec<u8> {
    if seqs.is_empty() {
        return Vec::new();
    }
    let width = seqs[0].len();
    let mut consensus = Vec::with_capacity(width);
    for col in 0..width {
        let mut counts: HashMap<u8, usize> = HashMap::new();
        for seq in seqs {
            if col < seq.len() {
                let base = seq[col];
                if base != b'-' && base != b'.' {
                    *counts.entry(base).or_insert(0) += 1;
                }
            }
        }
        let total: usize = counts.values().sum();
        let best_count = counts.values().copied().max().unwrap_or(0);
        let mut best = counts
            .into_iter()
            .filter_map(|(base, count)| (count == best_count).then_some(base));
        let best_char = best.next();
        if total > 0 && best.next().is_none() && (best_count as f64) / (total as f64) >= 0.7 {
            consensus.push(best_char.unwrap());
        } else {
            consensus.push(b'X');
        }
    }
    consensus
}

fn run_muscle(muscle_bin: &str, records: &[FastaRecord]) -> anyhow::Result<Vec<Vec<u8>>> {
    let mut child = Command::new(muscle_bin)
        .arg("-clwstrict")
        .stdin(Stdio::piped())
        .stdout(Stdio::piped())
        .stderr(Stdio::piped())
        .spawn()?;
    {
        let stdin = child
            .stdin
            .as_mut()
            .ok_or_else(|| anyhow::anyhow!("failed to open MUSCLE stdin"))?;
        for record in records {
            writeln!(stdin, ">{}\n{}", record.id, record.sequence)?;
        }
    }
    let output = child.wait_with_output()?;
    anyhow::ensure!(
        output.status.success(),
        "MUSCLE failed with {}: {}",
        output.status,
        String::from_utf8_lossy(&output.stderr)
    );
    parse_clustal(&String::from_utf8(output.stdout)?)
}

fn parse_clustal(text: &str) -> anyhow::Result<Vec<Vec<u8>>> {
    let mut order = Vec::new();
    let mut sequences: HashMap<String, Vec<u8>> = HashMap::new();
    for line in text.lines() {
        if line.trim().is_empty()
            || line.starts_with("CLUSTAL")
            || line.starts_with("MUSCLE")
            || line.starts_with(char::is_whitespace)
        {
            continue;
        }
        let mut fields = line.split_whitespace();
        let Some(id) = fields.next() else { continue };
        let Some(sequence) = fields.next() else {
            continue;
        };
        if !sequences.contains_key(id) {
            order.push(id.to_string());
        }
        sequences
            .entry(id.to_string())
            .or_default()
            .extend_from_slice(sequence.as_bytes());
    }
    anyhow::ensure!(!order.is_empty(), "MUSCLE returned no Clustal sequences");
    let aligned = order
        .iter()
        .map(|id| sequences.remove(id).unwrap_or_default())
        .collect::<Vec<_>>();
    let width = aligned[0].len();
    anyhow::ensure!(
        aligned.iter().all(|sequence| sequence.len() == width),
        "MUSCLE returned unequal Clustal sequence lengths"
    );
    Ok(aligned)
}

#[allow(clippy::too_many_arguments)]
pub fn run(
    input: &Path,
    output: &Path,
    muscle_bin: Option<&str>,
    mafft_bin: Option<&str>,
    per_locus_dir: Option<&Path>,
) -> anyhow::Result<()> {
    let records = phyluce_io::read_fasta(input)
        .with_context(|| format!("reading input FASTA {}", input.display()))?;
    crate::cli_warn!(
        "There are {} baits in {}",
        records.len(),
        input.file_name().and_then(|s| s.to_str()).unwrap_or("")
    );

    let mut buckets: Vec<String> = Vec::new();
    let mut d: HashMap<String, Vec<FastaRecord>> = HashMap::new();
    for record in records {
        let locus = record
            .id
            .split('_')
            .next()
            .unwrap_or(&record.id)
            .to_string();
        if !d.contains_key(&locus) {
            buckets.push(locus.clone());
        }
        d.entry(locus).or_default().push(record);
    }

    let mut out = std::fs::File::create(output)
        .with_context(|| format!("creating output file {}", output.display()))?;
    if let Some(dir) = per_locus_dir {
        std::fs::create_dir_all(dir)
            .with_context(|| format!("creating per-locus output directory {}", dir.display()))?;
    }
    let mut count = 0usize;
    for locus in &buckets {
        let recs = &d[locus];
        let sequence = if recs.len() > 1 {
            let seqs = if let Some(muscle_bin) = muscle_bin {
                run_muscle(muscle_bin, recs)?
            } else {
                let resolved;
                let mafft_bin = match mafft_bin {
                    Some(binary) => binary,
                    None => {
                        let config = phyluce_config::PhyluceConfig::load()?;
                        resolved = config.get_user_path("binaries", "mafft")?;
                        &resolved
                    }
                };
                let inputs: Vec<(String, String)> = recs
                    .iter()
                    .map(|r| (r.id.clone(), r.sequence.clone()))
                    .collect();
                let alignment = phyluce_align::mafft::run_mafft(mafft_bin, &inputs)?;
                alignment.rows.into_iter().map(|row| row.seq).collect()
            };
            String::from_utf8_lossy(&dumb_consensus(&seqs)).into_owned()
        } else {
            recs[0].sequence.clone()
        };
        writeln!(out, ">{locus}\n{sequence}")?;
        if let Some(dir) = per_locus_dir {
            // One unwrapped, single-record FASTA per locus, named to match
            // GeneMiner2-UCE's `-r` reference-directory convention
            // (`<locus>.fasta`, e.g. `uce-0.fasta`).
            let locus_path = dir.join(format!("{locus}.fasta"));
            std::fs::write(&locus_path, format!(">{locus}\n{sequence}\n"))
                .with_context(|| format!("writing per-locus FASTA {}", locus_path.display()))?;
        }
        count += 1;
    }

    crate::cli_warn!(
        "Wrote {count} loci to {}",
        output.file_name().and_then(|s| s.to_str()).unwrap_or("")
    );
    if let Some(dir) = per_locus_dir {
        crate::cli_warn!("Wrote {count} per-locus FASTA files to {}", dir.display());
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn per_locus_dir_writes_one_unwrapped_file_per_locus() {
        // Single-probe loci only, so this doesn't need a real MAFFT/MUSCLE
        // binary on PATH -- exercises the per-locus-dir plumbing (dir
        // creation, one file per locus, unwrapped single-line content,
        // `<locus>.fasta` naming) independent of the consensus step.
        let dir = std::env::temp_dir().join(format!(
            "phyluce-reconstruct-per-locus-{}",
            std::process::id()
        ));
        let _ = std::fs::remove_dir_all(&dir);
        std::fs::create_dir_all(&dir).unwrap();
        let input = dir.join("probes.fasta");
        std::fs::write(
            &input,
            ">uce-1_p1 |design:test\nACGTACGT\n>uce-2_p1 |design:test\nTTTTGGGG\n",
        )
        .unwrap();
        let output = dir.join("combined.fasta");
        let per_locus_dir = dir.join("per_locus");

        run(&input, &output, None, None, Some(&per_locus_dir)).unwrap();

        let uce1 = std::fs::read_to_string(per_locus_dir.join("uce-1.fasta")).unwrap();
        assert_eq!(uce1, ">uce-1\nACGTACGT\n");
        let uce2 = std::fs::read_to_string(per_locus_dir.join("uce-2.fasta")).unwrap();
        assert_eq!(uce2, ">uce-2\nTTTTGGGG\n");

        std::fs::remove_dir_all(&dir).unwrap();
    }

    #[test]
    fn consensus_matches_biopython_gap_and_threshold_rules() {
        assert_eq!(dumb_consensus(&[b"A-C".to_vec(), b"A.C".to_vec()]), b"AXC");
        assert_eq!(
            dumb_consensus(&[b"A".to_vec(), b"-".to_vec(), b".".to_vec()]),
            b"A"
        );
        assert_eq!(
            dumb_consensus(&[b"A".to_vec(), b"A".to_vec(), b"G".to_vec()]),
            b"X"
        );
    }

    #[test]
    fn parses_interleaved_clustal_output() {
        let text = "CLUSTAL W\n\na  AC-\nb  ACG\n   ** \n\na  GT\nb  GT\n";
        assert_eq!(parse_clustal(text).unwrap(), [b"AC-GT", b"ACGGT"]);
    }
}
