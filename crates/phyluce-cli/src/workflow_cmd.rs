//! CLI wiring for `phyluce workflow run`, mirroring `phyluce_workflow`.
//!
//! The Python original calls Snakemake's Python API directly
//! (`snakemake.snakemake(...)`); this shells out to the `snakemake` CLI
//! binary instead with the equivalent flags. Untested (Snakemake isn't
//! installed in this environment).

use std::path::Path;

use phyluce_config::PhyluceConfig;

pub fn run(
    config: &Path,
    output: &Path,
    workflow: &str,
    cores: u32,
    dryrun: bool,
) -> anyhow::Result<()> {
    anyhow::ensure!(
        matches!(workflow, "mapping" | "correction" | "phasing"),
        "--workflow must be one of: mapping, correction, phasing"
    );
    std::fs::create_dir_all(output)?;

    let cfg = PhyluceConfig::load()?;
    let snake_file = cfg.get_user_path("workflows", workflow)?;
    let snakemake_bin = cfg
        .get_user_path("binaries", "snakemake")
        .unwrap_or_else(|_| "snakemake".to_string());

    let mut cmd = std::process::Command::new(&snakemake_bin);
    cmd.arg("--snakefile")
        .arg(&snake_file)
        .arg("--cores")
        .arg(cores.to_string())
        .arg("--configfile")
        .arg(config)
        .arg("--directory")
        .arg(output);
    if dryrun {
        cmd.arg("--dry-run");
    }
    let status = cmd.status()?;
    anyhow::ensure!(status.success(), "snakemake exited with {status}");
    Ok(())
}
