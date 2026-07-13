//! CLI wiring for `phyluce genetrees sort-multilocus-bootstraps`, mirroring
//! `phyluce_genetrees_sort_multilocus_bootstraps`.
//!
//! Reads the plain-text replicate format written by
//! `phyluce genetrees generate-multilocus-bootstrap-count` (not Python
//! `pickle` -- see `phyluce_genetrees::bootstrap`'s docs).

use std::collections::HashMap;
use std::path::Path;

pub fn run(input: &Path, bootstrap_replicates: &Path, output: &Path) -> anyhow::Result<()> {
    let replicates = phyluce_genetrees::bootstrap::read_replicates(bootstrap_replicates)?;

    crate::cli_info!("Reading bootstrap replicates");
    let mut all_bootreps: HashMap<String, Vec<String>> = HashMap::new();
    for entry in std::fs::read_dir(input)? {
        let dir = entry?.path();
        if !dir.is_dir() {
            continue;
        }
        let dir_name = dir
            .file_name()
            .and_then(|n| n.to_str())
            .unwrap_or("")
            .to_string();
        let mut bootrep_files: Vec<std::path::PathBuf> = std::fs::read_dir(&dir)?
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
        let text = std::fs::read_to_string(bootrep_files.remove(0))?;
        let lines: Vec<String> = text.lines().map(|l| format!("{l}\n")).collect();
        all_bootreps.insert(dir_name, lines);
        print!(".");
    }

    crate::cli_info!("\nWriting bootstrap replicates");
    std::fs::create_dir_all(output)?;
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
            anyhow::ensure!(
                !lines.is_empty(),
                "ran out of bootstrap replicates for locus {locus_name}"
            );
            out.push_str(&lines.remove(0));
        }
        std::fs::write(out_path, out)?;
        print!(".");
    }
    crate::cli_info!();
    Ok(())
}
