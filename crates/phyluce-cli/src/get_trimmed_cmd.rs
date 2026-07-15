//! CLI wiring for `phyluce align get-trimmed-alignments-from-untrimmed`,
//! mirroring `phyluce_align_get_trimmed_alignments_from_untrimmed`.

use std::path::{Path, PathBuf};

use phyluce_align::trim::{trim_alignment_running, validate_trim_parameters};
use phyluce_align::{nexus::format_nexus, Alignment};
use phyluce_io::read_fasta;

const FASTA_EXTENSIONS: &[&str] = &[".fasta", ".fsa", ".aln", ".fa"];

#[allow(clippy::too_many_arguments)]
pub fn run(
    alignments_dir: &Path,
    output_dir: &Path,
    window: usize,
    proportion: f64,
    threshold: f64,
    max_divergence: f64,
    min_length: usize,
    cores: usize,
) -> anyhow::Result<()> {
    anyhow::ensure!(cores > 0, "--cores must be greater than zero");
    validate_trim_parameters(window, proportion, threshold, max_divergence)?;
    crate::output_path::prepare_output_dir(output_dir)?;

    let mut files: Vec<PathBuf> = Vec::new();
    for entry in std::fs::read_dir(alignments_dir)? {
        let path = entry?.path();
        if let Some(name) = path.file_name().and_then(|n| n.to_str()) {
            if FASTA_EXTENSIONS.iter().any(|ext| name.ends_with(ext)) {
                files.push(path);
            }
        }
    }
    files.sort();
    crate::parallel::ensure_unique_output_names(
        files
            .iter()
            .map(|file| format!("{}.nexus", strip_known_extension(file))),
    )?;

    let results = crate::parallel::try_map_ordered(files, cores, |file| {
        let stem = strip_known_extension(&file);
        let records = read_fasta(&file)?;
        let alignment =
            Alignment::from_pairs(records.into_iter().map(|r| (r.id, r.sequence)).collect());
        alignment.validate()?;
        let trimmed = trim_alignment_running(
            &alignment,
            window,
            proportion,
            threshold,
            max_divergence,
            min_length,
        );
        match trimmed {
            Some(aln) => {
                let out_path = output_dir.join(format!("{stem}.nexus"));
                std::fs::write(out_path, format_nexus(&aln))?;
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

fn strip_known_extension(path: &Path) -> String {
    let name = path.file_name().and_then(|n| n.to_str()).unwrap_or("");
    for ext in FASTA_EXTENSIONS {
        if let Some(stripped) = name.strip_suffix(ext) {
            return stripped.to_string();
        }
    }
    path.file_stem()
        .and_then(|s| s.to_str())
        .unwrap_or(name)
        .to_string()
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn rejects_unequal_fasta_rows_without_panicking() {
        let root = std::env::temp_dir().join(format!(
            "phyluce-trim-invalid-alignment-{}",
            std::process::id()
        ));
        let input = root.join("input");
        let output = root.join("output");
        let _ = std::fs::remove_dir_all(&root);
        std::fs::create_dir_all(&input).unwrap();
        std::fs::write(input.join("bad.fasta"), ">a\nAAAA\n>b\nAA\n").unwrap();

        let error = run(&input, &output, 5, 0.65, 0.65, 0.2, 1, 2).unwrap_err();
        assert!(error.to_string().contains("expected 4"));
        std::fs::remove_dir_all(root).unwrap();
    }
}
