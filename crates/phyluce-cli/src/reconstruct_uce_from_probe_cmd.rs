//! CLI wiring for `phyluce probe reconstruct-uce-from-probe`, mirroring
//! `phyluce_probe_reconstruct_uce_from_probe`.
//!
//! The Python original aligns multi-probe loci with MUSCLE + Clustal
//! format and takes Biopython's `AlignInfo.dumb_consensus()`. MUSCLE isn't
//! available in this environment, so this port aligns with MAFFT instead
//! (`phyluce_align::mafft::run_mafft`) and re-implements
//! `dumb_consensus`'s per-column majority-vote logic directly -- a
//! deliberate, documented divergence, not a bug. Single-probe loci (the
//! common case) need no external tool at all and match byte-for-byte.

use std::collections::HashMap;
use std::io::Write as _;
use std::path::Path;

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
                *counts.entry(seq[col].to_ascii_uppercase()).or_insert(0) += 1;
            }
        }
        let total: usize = counts.values().sum();
        let (best_char, best_count) = counts
            .into_iter()
            .max_by_key(|&(_, c)| c)
            .unwrap_or((b'X', 0));
        if total > 0 && (best_count as f64) / (total as f64) >= 0.7 {
            consensus.push(best_char);
        } else {
            consensus.push(b'X');
        }
    }
    consensus
}

pub fn run(input: &Path, output: &Path, mafft_bin: Option<&str>) -> anyhow::Result<()> {
    let records = phyluce_io::read_fasta(input)?;
    eprintln!(
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

    let mut out = std::fs::File::create(output)?;
    let mut count = 0usize;
    for locus in &buckets {
        let recs = &d[locus];
        if recs.len() > 1 {
            let mafft_bin = mafft_bin.ok_or_else(|| {
                anyhow::anyhow!("locus '{locus}' has {} probes and needs alignment, but no mafft binary was given", recs.len())
            })?;
            let inputs: Vec<(String, String)> = recs
                .iter()
                .map(|r| (r.id.clone(), r.sequence.clone()))
                .collect();
            let alignment = phyluce_align::mafft::run_mafft(mafft_bin, &inputs)?;
            let seqs: Vec<Vec<u8>> = alignment.rows.iter().map(|r| r.seq.clone()).collect();
            let consensus = dumb_consensus(&seqs);
            writeln!(out, ">{locus}\n{}", String::from_utf8_lossy(&consensus))?;
        } else {
            writeln!(out, ">{locus}\n{}", recs[0].sequence)?;
        }
        count += 1;
    }

    eprintln!(
        "Wrote {count} loci to {}",
        output.file_name().and_then(|s| s.to_str()).unwrap_or("")
    );
    Ok(())
}
