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
) -> anyhow::Result<()> {
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

    let mut dropped = 0usize;
    for file in &files {
        let stem = strip_known_extension(file);
        let records = read_fasta(file)?;
        let alignment =
            Alignment::from_pairs(records.into_iter().map(|r| (r.id, r.sequence)).collect());
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
                print!(".");
                let out_path = output_dir.join(format!("{stem}.nexus"));
                std::fs::write(out_path, format_nexus(&aln))?;
            }
            None => {
                print!("X");
                dropped += 1;
            }
        }
        use std::io::Write as _;
        std::io::stdout().flush().ok();
    }
    crate::cli_info!();
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
