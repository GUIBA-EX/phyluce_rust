//! CLI wiring for `phyluce align reduce-alignments-with-raxml`, mirroring
//! `phyluce_align_reduce_alignments_with_raxml`.

use std::path::Path;

use anyhow::Context;
use phyluce_config::PhyluceConfig;
use phyluce_external::ExternalCommand;

use crate::informative_sites_cmd::find_alignment_files;

pub fn run(alignments_dir: &Path, output_dir: &Path, input_format: &str) -> anyhow::Result<()> {
    crate::output_path::prepare_output_dir(output_dir)?;
    let cfg = PhyluceConfig::load()?;
    let raxml_bin = cfg.get_user_path("binaries", "raxmlHPC-SSE3")?;

    let files = find_alignment_files(alignments_dir, input_format)
        .with_context(|| format!("reading alignments directory {}", alignments_dir.display()))?;
    for file in &files {
        let old_name = file.file_name().and_then(|n| n.to_str()).unwrap_or("");
        let tdir =
            std::env::temp_dir().join(format!("phyluce-raxml-reduce-{}", std::process::id()));
        std::fs::create_dir_all(&tdir)
            .with_context(|| format!("creating temp directory {}", tdir.display()))?;

        // `tdir` must be removed whether RAxML succeeds or fails (a
        // malformed/degenerate locus can make RAxML exit non-zero), so run
        // the fallible steps in a closure and always clean up afterward
        // rather than an early `?` return skipping the `remove_dir_all`.
        let result = (|| -> anyhow::Result<()> {
            ExternalCommand::new(&raxml_bin)
                .args([
                    "-f".to_string(),
                    "c".to_string(),
                    "-m".to_string(),
                    "GTRGAMMA".to_string(),
                    "-s".to_string(),
                    file.to_string_lossy().to_string(),
                    "-w".to_string(),
                    tdir.to_string_lossy().to_string(),
                    "-n".to_string(),
                    "test".to_string(),
                ])
                .run()
                .with_context(|| format!("running raxml for {}", file.display()))?;

            let reduced_pth = file.with_file_name(format!("{old_name}.reduced"));
            let new_pth = output_dir.join(old_name);
            if reduced_pth.exists() {
                std::fs::rename(&reduced_pth, &new_pth).with_context(|| {
                    format!(
                        "moving reduced alignment {} to {}",
                        reduced_pth.display(),
                        new_pth.display()
                    )
                })?;
            } else {
                std::fs::copy(file, &new_pth).with_context(|| {
                    format!("copying {} to {}", file.display(), new_pth.display())
                })?;
            }
            Ok(())
        })();
        let _ = std::fs::remove_dir_all(&tdir);
        result?;
        print!(".");
    }
    crate::cli_info!();
    Ok(())
}
