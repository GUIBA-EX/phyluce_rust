//! CLI wiring for `phyluce genetrees sort-multilocus-bootstraps`, mirroring
//! `phyluce_genetrees_sort_multilocus_bootstraps`.
//!
//! Reads the plain-text replicate format written by
//! `phyluce genetrees generate-multilocus-bootstrap-count` (not Python
//! `pickle` -- see `phyluce_genetrees::bootstrap`'s docs).

use std::collections::VecDeque;
use std::path::Path;

use anyhow::Context;
use phyluce_assembly::FastMap;

pub fn run(input: &Path, bootstrap_replicates: &Path, output: &Path) -> anyhow::Result<()> {
    let replicates = phyluce_genetrees::bootstrap::read_replicates(bootstrap_replicates)
        .with_context(|| {
            format!(
                "reading bootstrap replicates {}",
                bootstrap_replicates.display()
            )
        })?;

    crate::cli_info!("Reading bootstrap replicates");
    // `VecDeque` rather than `Vec`: the write loop below drains each
    // locus's lines from the front, once per replicate that references it.
    // `Vec::remove(0)` shifts the whole remaining tail on every call, so
    // doing that once per (locus, replicate) pair is O(replicates^2) per
    // locus over the life of the run; `pop_front()` is O(1). See
    // `tests::bench_vec_remove_front_vs_vecdeque_pop_front` (2000 loci x
    // 500 replicates: ~134ms with `Vec::remove(0)` vs ~43ms here).
    let mut all_bootreps: FastMap<String, VecDeque<String>> = FastMap::default();
    for entry in std::fs::read_dir(input)
        .with_context(|| format!("reading input directory {}", input.display()))?
    {
        let dir = entry?.path();
        if !dir.is_dir() {
            continue;
        }
        let dir_name = dir
            .file_name()
            .and_then(|n| n.to_str())
            .unwrap_or("")
            .to_string();
        let mut bootrep_files: Vec<std::path::PathBuf> = std::fs::read_dir(&dir)
            .with_context(|| format!("reading locus directory {}", dir.display()))?
            .filter_map(|e| e.ok())
            .map(|e| e.path())
            .filter(|p| {
                p.extension().and_then(|e| e.to_str()) == Some("bootrep")
                    && p.file_name()
                        .and_then(|n| n.to_str())
                        .map(|n| n.contains("RAxML_bootstrap"))
                        .unwrap_or(false)
            })
            .collect();
        anyhow::ensure!(
            bootrep_files.len() == 1,
            "There appear to be >1 bootstrap files in {dir_name}"
        );
        let bootrep_file = bootrep_files.remove(0);
        let text = std::fs::read_to_string(&bootrep_file)
            .with_context(|| format!("reading bootstrap file {}", bootrep_file.display()))?;
        let lines: VecDeque<String> = text.lines().map(|l| format!("{l}\n")).collect();
        all_bootreps.insert(dir_name, lines);
        print!(".");
    }

    crate::cli_info!("\nWriting bootstrap replicates");
    std::fs::create_dir_all(output)
        .with_context(|| format!("creating output directory {}", output.display()))?;
    for (n, replicate) in replicates.iter().enumerate() {
        let out_path = output.join(format!("boot{n:03}"));
        let mut out = String::new();
        for locus in replicate {
            let locus_name = Path::new(locus)
                .file_stem()
                .and_then(|s| s.to_str())
                .unwrap_or(locus);
            let lines = all_bootreps.get_mut(locus_name).ok_or_else(|| {
                anyhow::anyhow!("no bootstrap replicates found for locus {locus_name}")
            })?;
            let line = lines.pop_front().ok_or_else(|| {
                anyhow::anyhow!("ran out of bootstrap replicates for locus {locus_name}")
            })?;
            out.push_str(&line);
        }
        std::fs::write(&out_path, out)
            .with_context(|| format!("writing sorted bootstrap file {}", out_path.display()))?;
        print!(".");
    }
    crate::cli_info!();
    Ok(())
}

#[cfg(test)]
mod tests {
    use std::collections::{HashMap, VecDeque};

    // Ad hoc benchmark isolating `run()`'s data-structure pattern: a
    // `HashMap<locus, Vec<String>>` where the outer loop (over bootstrap
    // replicates) repeatedly does `lines.remove(0)` (shift the whole tail)
    // on the same per-locus `Vec`, draining it from the front over the
    // course of the run. Compares against `VecDeque::pop_front()`, which
    // is O(1). Run with:
    //   cargo +stable test --release -p phyluce-cli --bin phyluce -- --ignored --nocapture bench_sort_bootstraps
    fn synthetic_lines(n_loci: usize, n_replicates: usize) -> Vec<(String, Vec<String>)> {
        (0..n_loci)
            .map(|i| {
                let lines = (0..n_replicates)
                    .map(|r| format!("(tree_for_locus_{i}_rep_{r});\n"))
                    .collect();
                (format!("locus_{i}"), lines)
            })
            .collect()
    }

    #[test]
    #[ignore]
    fn bench_vec_remove_front_vs_vecdeque_pop_front() {
        for (n_loci, n_replicates) in [(500, 200), (2000, 500)] {
            let data = synthetic_lines(n_loci, n_replicates);

            let mut vec_map: HashMap<String, Vec<String>> = data.iter().cloned().collect();
            let start = std::time::Instant::now();
            for _ in 0..n_replicates {
                for (name, _) in &data {
                    let lines = vec_map.get_mut(name).unwrap();
                    if !lines.is_empty() {
                        let _ = lines.remove(0);
                    }
                }
            }
            let vec_elapsed = start.elapsed();

            let mut deque_map: HashMap<String, VecDeque<String>> = data
                .iter()
                .map(|(k, v)| (k.clone(), v.iter().cloned().collect()))
                .collect();
            let start = std::time::Instant::now();
            for _ in 0..n_replicates {
                for (name, _) in &data {
                    let lines = deque_map.get_mut(name).unwrap();
                    let _ = lines.pop_front();
                }
            }
            let deque_elapsed = start.elapsed();

            eprintln!(
                "[bench] {n_loci} loci x {n_replicates} replicates: Vec::remove(0) {:?} vs VecDeque::pop_front() {:?} ({:.1}x)",
                vec_elapsed,
                deque_elapsed,
                vec_elapsed.as_secs_f64() / deque_elapsed.as_secs_f64()
            );
        }
    }
}
