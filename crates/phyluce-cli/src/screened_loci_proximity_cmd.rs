//! CLI wiring for `phyluce probe get-screened-loci-by-proximity`, mirroring
//! `phyluce_probe_get_screened_loci_by_proximity`.
//!
//! The Python original breaks proximity-cluster ties with
//! `random.choice`, so it's inherently non-deterministic; this port
//! always keeps the *first* (lowest locus id) member of a cluster instead
//! of a random one -- a deliberate, documented divergence, not a bug.

use std::collections::{HashMap, HashSet};
use std::path::Path;

use phyluce_io::{read_fasta, write_fasta_record};
use regex::Regex;

fn parse_probe_position(locus: i64, field: &str) -> Option<(String, i64, i64)> {
    let mut chromo = None;
    let mut start = None;
    let mut end = None;
    for kv in field.split(',') {
        if let Some((k, v)) = kv.split_once(':') {
            match k {
                "probes-global-chromo" => chromo = Some(v.to_string()),
                "probes-global-start" => start = v.parse::<i64>().ok(),
                "probes-global-end" => end = v.parse::<i64>().ok(),
                _ => {}
            }
        }
    }
    let _ = locus;
    Some((chromo?, start?, end?))
}

/// Simple gap-based interval clustering (sort by start, merge intervals
/// whose gap to the running cluster is <= `distance`), standing in for
/// `bx.intervals.cluster.ClusterTree(distance, 2)`.
fn cluster_loci(positions: &[(i64, String, i64, i64)], distance: i64) -> Vec<HashSet<i64>> {
    let mut by_chromo: HashMap<&str, Vec<&(i64, String, i64, i64)>> = HashMap::new();
    for p in positions {
        by_chromo.entry(p.1.as_str()).or_default().push(p);
    }
    let mut clusters = Vec::new();
    for (_, mut items) in by_chromo {
        items.sort_by_key(|p| p.2);
        let mut current: Vec<&(i64, String, i64, i64)> = Vec::new();
        let mut current_end = i64::MIN;
        for item in items {
            if current.is_empty() || item.2 <= current_end + distance {
                current.push(item);
                current_end = current_end.max(item.3);
            } else {
                if current.len() > 1 {
                    clusters.push(current.iter().map(|p| p.0).collect());
                }
                current = vec![item];
                current_end = item.3;
            }
        }
        if current.len() > 1 {
            clusters.push(current.iter().map(|p| p.0).collect());
        }
    }
    clusters
}

pub fn run(input: &Path, output: &Path, distance: i64, regex_str: &str) -> anyhow::Result<()> {
    let regex = Regex::new(regex_str)?;
    let records = read_fasta(input)?;

    let mut positions_set: HashSet<(i64, String, i64, i64)> = HashSet::new();
    let mut starting_baits = 0usize;
    for record in &records {
        let locus: i64 = regex
            .captures(&record.id)
            .and_then(|c| c.get(1))
            .and_then(|m| m.as_str().parse().ok())
            .ok_or_else(|| anyhow::anyhow!("no regex match for probe id {:?}", record.id))?;
        let field =
            record.description.split('|').nth(1).ok_or_else(|| {
                anyhow::anyhow!("record '{}': missing '|' metadata field", record.id)
            })?;
        if let Some((chromo, start, end)) = parse_probe_position(locus, field) {
            positions_set.insert((locus, chromo, start, end));
        }
        starting_baits += 1;
    }
    let loci_count = positions_set
        .iter()
        .map(|p| p.0)
        .collect::<HashSet<_>>()
        .len();
    crate::cli_warn!("Start with {loci_count} loci and {starting_baits} baits");

    let positions: Vec<(i64, String, i64, i64)> = positions_set.into_iter().collect();
    let clusters = cluster_loci(&positions, distance);

    let mut bad: HashSet<i64> = HashSet::new();
    for cluster in &clusters {
        let mut loci: Vec<i64> = cluster.iter().cloned().collect();
        loci.sort();
        if loci.len() > 1 {
            // keep the first (lowest id); drop the rest -- see module docs
            // re: this being a deterministic stand-in for `random.choice`.
            bad.extend(loci.into_iter().skip(1));
        }
    }
    crate::cli_warn!(
        "Removing {} loci that appear to be within {distance} bp of one another",
        bad.len()
    );

    let mut out = std::fs::File::create(output)?;
    let mut kept_baits = 0usize;
    let mut kept_loci = HashSet::new();
    for record in &records {
        let locus: i64 = regex
            .captures(&record.id)
            .and_then(|c| c.get(1))
            .and_then(|m| m.as_str().parse().ok())
            .unwrap();
        if !bad.contains(&locus) {
            write_fasta_record(&mut out, &record.description, &record.sequence)?;
            kept_baits += 1;
            kept_loci.insert(locus);
        }
    }
    crate::cli_warn!("Ends with {} loci and {kept_baits} baits", kept_loci.len());
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn clusters_nearby_loci_and_leaves_distant_ones_alone() {
        let positions = vec![
            (1, "chr1".to_string(), 0, 100),
            (2, "chr1".to_string(), 150, 250), // within 10000bp of locus 1
            (3, "chr1".to_string(), 50000, 50100), // far away
        ];
        let clusters = cluster_loci(&positions, 10000);
        assert_eq!(clusters.len(), 1);
        assert_eq!(clusters[0], HashSet::from([1, 2]));
    }
}
