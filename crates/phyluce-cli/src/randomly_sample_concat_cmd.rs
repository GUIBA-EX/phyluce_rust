//! CLI wiring for `phyluce align randomly-sample-and-concatenate`,
//! mirroring `phyluce_align_randomly_sample_and_concatenate`.
//!
//! Uses a seeded xorshift PRNG for the without-replacement sample instead
//! of `numpy.random.choice`, consistent with the pattern in
//! `bootstrap_count_cmd.rs` -- reproducibility of Python's exact RNG
//! stream isn't the point, just an unbiased random sample each run.

use std::path::{Path, PathBuf};

use phyluce_align::concat::{concatenate, format_sets_block};
use phyluce_align::nexus::format_nexus_with_interleave;

use crate::informative_sites_cmd::load_alignment;

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

    /// Fisher-Yates partial shuffle, taking the first `size` elements --
    /// an unbiased without-replacement sample.
    fn sample_without_replacement(&mut self, items: &mut [PathBuf], size: usize) {
        let n = items.len();
        for i in 0..size.min(n) {
            let j = i + (self.next_u64() as usize) % (n - i);
            items.swap(i, j);
        }
    }
}

pub fn run(
    alignments_dir: &Path,
    output_dir: &Path,
    sample_size: usize,
    replicates: usize,
) -> anyhow::Result<()> {
    crate::output_path::prepare_output_dir(output_dir)?;
    let mut files: Vec<PathBuf> = std::fs::read_dir(alignments_dir)?
        .filter_map(|e| e.ok())
        .map(|e| e.path())
        .filter(|p| {
            p.file_name()
                .and_then(|n| n.to_str())
                .is_some_and(|n| n.contains(".nex"))
        })
        .collect();
    anyhow::ensure!(
        files.len() >= sample_size,
        "You have requested a larger sample size than your population of alignments"
    );

    let mut rng = SimpleRng::new();
    for i in 0..replicates {
        rng.sample_without_replacement(&mut files, sample_size);
        let sample = &files[..sample_size];

        let mut loaded: Vec<(String, phyluce_align::Alignment)> = Vec::with_capacity(sample.len());
        for file in sample {
            let name = file
                .file_name()
                .and_then(|n| n.to_str())
                .unwrap_or("")
                .to_string();
            let alignment = load_alignment(file, "nexus")?;
            loaded.push((name, alignment));
        }
        let (combined, charsets) = concatenate(&loaded);
        let mut text = format_nexus_with_interleave(&combined, false);
        text.push_str(&format_sets_block(&charsets));

        let align_name = format!("random-sample-{i}_{sample_size}-loci.nex");
        std::fs::write(output_dir.join(align_name), text)?;
    }
    Ok(())
}
