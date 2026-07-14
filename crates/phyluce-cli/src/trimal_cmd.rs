//! CLI wiring for `phyluce align get-trimal-trimmed-alignments-from-untrimmed`,
//! mirroring `phyluce_align_get_trimal_trimmed_alignments_from_untrimmed`.
//!
//! Untested against a live trimAl binary in this environment (not
//! installed) -- ported carefully from source, treat as best-effort until
//! validated against a real run.

use std::path::Path;

use anyhow::Context;
use phyluce_align::nexus::format_nexus;
use phyluce_config::PhyluceConfig;
use phyluce_io::read_fasta;

use crate::informative_sites_cmd::find_alignment_files;

pub fn run(
    alignments_dir: &Path,
    output_dir: &Path,
    input_format: &str,
    output_format: &str,
) -> anyhow::Result<()> {
    anyhow::ensure!(
        matches!(output_format, "fasta" | "nexus"),
        "output format '{output_format}' is not supported (only fasta/nexus)"
    );
    crate::output_path::prepare_output_dir(output_dir)?;

    let cfg = PhyluceConfig::load()?;
    let trimal_bin = cfg.get_user_path("binaries", "trimal")?;

    let files = find_alignment_files(alignments_dir, input_format)?;
    for file in &files {
        let name = file
            .file_name()
            .and_then(|n| n.to_str())
            .unwrap_or("")
            .split('.')
            .next()
            .unwrap_or("")
            .to_string();

        let trimmed_path = std::path::PathBuf::from(format!("{}-trimal", file.display()));
        let output = std::process::Command::new(&trimal_bin)
            .arg("-in")
            .arg(file)
            .arg("-out")
            .arg(&trimmed_path)
            .arg("-automated1")
            .arg("-fasta")
            .output()
            .with_context(|| format!("running trimAl for locus {name}"))?;
        anyhow::ensure!(
            output.status.success(),
            "trimAl failed for locus {name}: {}",
            String::from_utf8_lossy(&output.stderr).trim()
        );

        if !trimmed_path.is_file() {
            anyhow::bail!("trimAl did not create output for locus {name}");
        }
        let records = read_fasta(&trimmed_path)?;
        std::fs::remove_file(&trimmed_path)?;
        if records.is_empty() {
            crate::cli_warn!("Missing information for locus {name}");
            print!(".");
            continue;
        }
        let trimmed = phyluce_align::Alignment::from_pairs(
            records.into_iter().map(|r| (r.id, r.sequence)).collect(),
        );

        let ext = output_format;
        let out_path = output_dir.join(format!("{name}.{ext}"));
        if output_format == "fasta" {
            let mut out = std::fs::File::create(out_path)?;
            for row in &trimmed.rows {
                phyluce_io::write_fasta_record(&mut out, &row.id, std::str::from_utf8(&row.seq)?)?;
            }
        } else {
            std::fs::write(out_path, format_nexus(&trimmed))?;
        }
        print!(".");
    }
    crate::cli_info!();
    Ok(())
}
