//! CLI wiring for `phyluce assembly assemblo-spades`, mirroring
//! `phyluce_assembly_assemblo_spades`.
//!
//! Untested against a live SPAdes binary in this environment (not
//! installed) -- ported carefully from source, treat as best-effort until
//! validated against a real run.

use std::path::{Path, PathBuf};

use anyhow::Context;
use phyluce_assembly::raw_reads::get_input_files;
use phyluce_config::PhyluceConfig;

pub(crate) fn get_input_data_pub(
    config: Option<&Path>,
    dir: Option<&Path>,
) -> anyhow::Result<Vec<(String, PathBuf)>> {
    get_input_data(config, dir)
}

fn get_input_data(
    config: Option<&Path>,
    dir: Option<&Path>,
) -> anyhow::Result<Vec<(String, PathBuf)>> {
    if let Some(config) = config {
        let text = std::fs::read_to_string(config)?;
        let sections = crate::conf::parse_ini(&text);
        let entries = sections
            .get("samples")
            .ok_or_else(|| anyhow::anyhow!("no [samples] section in --config"))?;
        entries
            .iter()
            .map(|(name, path)| {
                let expanded = shellexpand_home(path);
                anyhow::ensure!(
                    Path::new(&expanded).is_dir(),
                    "{expanded} is not a directory"
                );
                Ok((name.clone(), PathBuf::from(expanded)))
            })
            .collect()
    } else if let Some(dir) = dir {
        let mut groups = Vec::new();
        for entry in std::fs::read_dir(dir)? {
            let path = entry?.path();
            let name = path.file_name().unwrap().to_string_lossy().to_string();
            groups.push((name, path));
        }
        groups.sort();
        Ok(groups)
    } else {
        anyhow::bail!("one of --config or --dir is required")
    }
}

fn shellexpand_home(path: &str) -> String {
    if let Some(rest) = path.strip_prefix('~') {
        if let Ok(home) = std::env::var("HOME") {
            return format!("{home}{rest}");
        }
    }
    path.to_string()
}

fn cleanup_assembly_directory(dir: &Path) -> anyhow::Result<()> {
    let names: Vec<String> = std::fs::read_dir(dir)?
        .filter_map(|e| e.ok())
        .map(|e| e.file_name().to_string_lossy().to_string())
        .collect();
    if !names.contains(&"scaffolds.fasta".to_string())
        || !names.contains(&"contigs.fasta".to_string())
        || !names.contains(&"spades.log".to_string())
    {
        crate::cli_warn!("Expected assembly files were not found in output.");
        return Ok(());
    }
    for entry in std::fs::read_dir(dir)? {
        let path = entry?.path();
        let name = path.file_name().unwrap().to_string_lossy().to_string();
        if !matches!(
            name.as_str(),
            "scaffolds.fasta" | "contigs.fasta" | "spades.log"
        ) {
            if path.is_dir() {
                std::fs::remove_dir_all(&path)?;
            } else {
                std::fs::remove_file(&path)?;
            }
        }
    }
    Ok(())
}

fn generate_symlink(contig_dir: &Path, sample_dir: &Path, sample: &str) -> anyhow::Result<()> {
    let spades_fname = sample_dir.join("contigs.fasta");
    anyhow::ensure!(
        spades_fname.is_file(),
        "SPAdes did not produce {}",
        spades_fname.display()
    );
    let relpth = pathdiff(&spades_fname, contig_dir)
        .with_context(|| format!("computing contig link for {sample}"))?;
    let link_path =
        crate::output_path::output_file(contig_dir, &format!("{sample}.contigs.fasta"))?;
    #[cfg(unix)]
    std::os::unix::fs::symlink(relpth, link_path)
        .with_context(|| format!("creating contig link for {sample}"))?;
    #[cfg(not(unix))]
    anyhow::bail!("contig symlinks are only supported on Unix");
    Ok(())
}

/// Minimal relative-path computation (std lacks one): assumes both paths
/// are absolute and share some common ancestor.
pub(crate) fn pathdiff_pub(target: &Path, base: &Path) -> Option<PathBuf> {
    pathdiff(target, base).ok()
}

fn pathdiff(target: &Path, base: &Path) -> std::io::Result<PathBuf> {
    let target_components: Vec<_> = target.components().collect();
    let base_components: Vec<_> = base.components().collect();
    let common = target_components
        .iter()
        .zip(base_components.iter())
        .take_while(|(a, b)| a == b)
        .count();
    let mut result = PathBuf::new();
    for _ in common..base_components.len() {
        result.push("..");
    }
    for c in &target_components[common..] {
        result.push(c.as_os_str());
    }
    Ok(result)
}

#[allow(clippy::too_many_arguments)]
pub fn run(
    output: &Path,
    cores: u32,
    memory: u32,
    subfolder: &str,
    no_clean: bool,
    config: Option<&Path>,
    dir: Option<&Path>,
) -> anyhow::Result<()> {
    let cfg = PhyluceConfig::load()?;
    let spades_bin = cfg.get_user_path("binaries", "spades").map_err(|_| {
        anyhow::anyhow!(
            "Cannot find SPAdes. Ensure $PATH is correctly entered in your phyluce.conf file."
        )
    })?;
    let cov_cutoff = cfg.get_user_param("spades", "cov_cutoff").unwrap_or("5");

    std::fs::create_dir_all(output)?;
    let contig_dir = output.join("contigs");
    std::fs::create_dir_all(&contig_dir)?;

    let input_data = get_input_data(config, dir)?;
    for (sample, sample_input_dir) in &input_data {
        crate::cli_info!("Processing {sample}");
        let sample_dir = crate::output_path::output_file(output, &format!("{sample}_spades"))?;
        std::fs::create_dir_all(&sample_dir)?;

        let reads = get_input_files(sample_input_dir, subfolder)?;
        match (&reads.r1, &reads.r2) {
            (Some(r1), Some(r2)) => {
                crate::cli_info!("Running SPAdes for PE data");
                let mut cmd = std::process::Command::new(&spades_bin);
                cmd.arg("--careful")
                    .arg("--sc")
                    .arg("--memory")
                    .arg(memory.to_string())
                    .arg("--threads")
                    .arg(cores.to_string())
                    .arg("--cov-cutoff")
                    .arg(cov_cutoff)
                    .arg("--pe1-1")
                    .arg(r1)
                    .arg("--pe1-2")
                    .arg(r2)
                    .arg("-o")
                    .arg(&sample_dir);
                if let Some(s) = &reads.singleton {
                    cmd.arg("--pe1-s").arg(s);
                }
                let status = cmd
                    .status()
                    .with_context(|| format!("running SPAdes for sample {sample}"))?;
                anyhow::ensure!(
                    status.success(),
                    "SPAdes failed for sample {sample}: {status}"
                );

                if !no_clean {
                    cleanup_assembly_directory(&sample_dir)?;
                }
                generate_symlink(&contig_dir, &sample_dir, sample)?;
            }
            (Some(_), None) => {
                crate::cli_warn!("assemblo-spades will not run single-end data");
            }
            _ => {}
        }
    }
    Ok(())
}
