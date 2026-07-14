//! CLI wiring for `phyluce align reduce-alignments-with-raxml`, mirroring
//! `phyluce_align_reduce_alignments_with_raxml`.
//!
//! Untested: `raxmlHPC-SSE3` isn't installed in this environment.

use std::path::Path;

use phyluce_config::PhyluceConfig;
use phyluce_external::ExternalCommand;

use crate::informative_sites_cmd::find_alignment_files;

pub fn run(alignments_dir: &Path, output_dir: &Path, input_format: &str) -> anyhow::Result<()> {
    crate::output_path::prepare_output_dir(output_dir)?;
    let cfg = PhyluceConfig::load()?;
    let raxml_bin = cfg.get_user_path("binaries", "raxmlHPC-SSE3")?;

    let files = find_alignment_files(alignments_dir, input_format)?;
    for file in &files {
        let old_name = file.file_name().and_then(|n| n.to_str()).unwrap_or("");
        let tdir =
            std::env::temp_dir().join(format!("phyluce-raxml-reduce-{}", std::process::id()));
        std::fs::create_dir_all(&tdir)?;

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
            .run()?;

        let reduced_pth = file.with_file_name(format!("{old_name}.reduced"));
        let new_pth = output_dir.join(old_name);
        if reduced_pth.exists() {
            std::fs::rename(&reduced_pth, &new_pth)?;
        } else {
            std::fs::copy(file, &new_pth)?;
        }
        std::fs::remove_dir_all(&tdir)?;
        print!(".");
    }
    crate::cli_info!();
    Ok(())
}
