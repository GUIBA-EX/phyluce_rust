//! CLI wiring for `phyluce genetrees generate-multilocus-bootstrap-count`,
//! mirroring `phyluce_genetrees_generate_multilocus_bootstrap_count`.
//!
//! Uses `phyluce_genetrees::bootstrap`'s plain-text replicate format
//! instead of Python `pickle` -- see that module's docs.

use std::collections::HashMap;
use std::path::Path;

const ALIGNMENT_EXTENSIONS: &[&str] = &[".phylip", ".phy", ".phylip-relaxed"];

/// A small xorshift-style PRNG, seeded from the OS, standing in for
/// Python's `random.choice` (exact reproducibility of Python's Mersenne
/// Twister stream isn't the point here -- the two ends of this pipeline
/// only need to agree with *each other*, and both ends are this same
/// Rust implementation; see `phyluce-genetrees::bootstrap`'s docs).
struct SimpleRng(u64);

impl SimpleRng {
    fn new() -> Self {
        let seed = std::time::SystemTime::now()
            .duration_since(std::time::UNIX_EPOCH)
            .map(|d| d.as_nanos() as u64)
            .unwrap_or(0x9E3779B97F4A7C15)
            | 1;
        SimpleRng(seed)
    }

    fn next_u64(&mut self) -> u64 {
        let mut x = self.0;
        x ^= x << 13;
        x ^= x >> 7;
        x ^= x << 17;
        self.0 = x;
        x
    }

    fn choice<'a>(&mut self, items: &'a [String]) -> &'a String {
        &items[(self.next_u64() as usize) % items.len()]
    }
}

pub fn run(
    alignments_dir: &Path,
    bootstrap_replicates: &Path,
    directory: &str,
    bootstrap_counts: &Path,
    bootreps: usize,
) -> anyhow::Result<()> {
    let mut loci: Vec<String> = Vec::new();
    for entry in std::fs::read_dir(alignments_dir)? {
        let path = entry?.path();
        if let Some(name) = path.file_name().and_then(|n| n.to_str()) {
            if ALIGNMENT_EXTENSIONS.iter().any(|ext| name.ends_with(ext)) {
                loci.push(name.to_string());
            }
        }
    }
    loci.sort();
    crate::cli_info!("Processing {} alignments", loci.len());

    let mut rng = SimpleRng::new();
    let replicates: Vec<Vec<String>> = (0..bootreps)
        .map(|_| (0..loci.len()).map(|_| rng.choice(&loci).clone()).collect())
        .collect();

    let mut counter: HashMap<&str, usize> = HashMap::new();
    for replicate in &replicates {
        for locus in replicate {
            *counter.entry(locus.as_str()).or_insert(0) += 1;
        }
    }

    let path = if directory.is_empty() {
        alignments_dir
    } else {
        Path::new(directory)
    };
    let mut out = String::new();
    for (locus, count) in &counter {
        out.push_str(&format!("{} {count}\n", path.join(locus).display()));
    }
    std::fs::write(bootstrap_counts, out)?;

    phyluce_genetrees::bootstrap::write_replicates(bootstrap_replicates, &replicates)?;
    Ok(())
}
