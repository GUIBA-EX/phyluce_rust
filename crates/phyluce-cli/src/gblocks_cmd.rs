//! CLI wiring for `phyluce align get-gblocks-trimmed-alignments-from-untrimmed`,
//! mirroring `phyluce_align_get_gblocks_trimmed_alignments_from_untrimmed`.
//!
//! Untested against a live Gblocks binary in this environment (not
//! installed/available on Apple Silicon per the legacy script's own
//! platform check) -- ported carefully from source, but treat as
//! best-effort until validated against a real run.

use std::path::Path;

use anyhow::Context;
use phyluce_align::nexus::format_nexus;
use phyluce_config::PhyluceConfig;
use phyluce_io::read_fasta;

use crate::informative_sites_cmd::{find_alignment_files, load_alignment};

#[allow(clippy::too_many_arguments)]
pub fn run(
    alignments_dir: &Path,
    output_dir: &Path,
    input_format: &str,
    b1: f64,
    b2: f64,
    b3: u32,
    b4: u32,
    output_format: &str,
) -> anyhow::Result<()> {
    anyhow::ensure!(b4 >= 2, "--b4 must be >= 2");
    std::fs::create_dir_all(output_dir)?;

    let cfg = PhyluceConfig::load()?;
    let gblocks_bin = cfg.get_user_path("binaries", "gblocks")?;

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
        let alignment = load_alignment(file, input_format)?;
        let taxa = alignment.ntax();

        let b1_arg = (b1 * taxa as f64).round() as i64 + 1;
        let mut b2_arg = (b2 * taxa as f64).round() as i64;
        if b2_arg < b1_arg {
            b2_arg = b1_arg;
        }

        let output = std::process::Command::new(&gblocks_bin)
            .arg(file)
            .arg("-t=DNA")
            .arg(format!("-b1={b1_arg}"))
            .arg(format!("-b2={b2_arg}"))
            .arg(format!("-b3={b3}"))
            .arg(format!("-b4={b4}"))
            .arg("-b5=h")
            .arg("-p=n")
            .output()
            .with_context(|| format!("running Gblocks for locus {name}"))?;
        // Gblocks conventionally exits non-zero even on success; the
        // legacy script never checks the exit code, only whether the
        // `-gb` output file exists afterward.

        // Mirrors `"{}-gb".format(align_file)`: a literal `-gb` suffix
        // appended to the whole path, not an extension replacement.
        let trimmed_path = std::path::PathBuf::from(format!("{}-gb", file.display()));
        if !trimmed_path.is_file() {
            anyhow::ensure!(
                output.status.success(),
                "Gblocks failed for locus {name}: {}",
                String::from_utf8_lossy(&output.stderr).trim()
            );
            crate::cli_warn!("Missing information for locus {name}");
            print!(".");
            continue;
        }
        let records = read_fasta(&trimmed_path)?;
        let trimmed = phyluce_align::Alignment::from_pairs(
            records.into_iter().map(|r| (r.id, r.sequence)).collect(),
        );
        std::fs::remove_file(&trimmed_path)?;

        let ext = if output_format == "fasta" {
            "fasta"
        } else {
            "nexus"
        };
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
