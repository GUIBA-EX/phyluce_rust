//! CLI wiring for `phyluce assembly assemblo-velvet`, mirroring
//! `phyluce_assembly_assemblo_velvet`.
//!
//! Untested against live velveth/velvetg binaries in this environment
//! (not installed) -- ported carefully from source, treat as best-effort
//! until validated against a real run.

use std::path::{Path, PathBuf};

use anyhow::Context;
use phyluce_assembly::raw_reads::{get_input_files, FileKind, ReadSet};
use phyluce_config::PhyluceConfig;

use crate::assemblo_spades_cmd::get_input_data_pub as get_input_data;

fn velvet_read_flag(reads: &ReadSet) -> &'static str {
    match (&reads.kind, reads.gzip) {
        (Some(FileKind::Fastq), true) => "-fastq.gz",
        (Some(FileKind::Fastq), false) => "-fastq",
        (Some(FileKind::Fasta), true) => "-fasta.gz",
        (Some(FileKind::Fasta), false) => "-fasta",
        (None, _) => "-fastq",
    }
}

fn run_velveth(velveth: &str, kmer: u32, reads: &ReadSet, output: &Path) -> anyhow::Result<()> {
    let name = format!("out_k{kmer}");
    let flag = velvet_read_flag(reads);
    let mut cmd = std::process::Command::new(velveth);
    cmd.current_dir(output)
        .arg(&name)
        .arg(kmer.to_string())
        .arg(flag)
        .arg("-separate")
        .arg("-shortPaired")
        .arg(reads.r1.as_ref().unwrap())
        .arg(reads.r2.as_ref().unwrap());
    if let Some(s) = &reads.singleton {
        cmd.arg("-short").arg(s);
    }
    let out = std::fs::File::create(output.join(format!("velveth-k{kmer}.out.log")))?;
    let err = std::fs::File::create(output.join(format!("velveth-k{kmer}.err.log")))?;
    cmd.stdout(out).stderr(err);
    let status = cmd.status().context("running velveth")?;
    anyhow::ensure!(status.success(), "velveth failed: {status}");
    Ok(())
}

fn run_velvetg(velvetg: &str, kmer: u32, output: &Path) -> anyhow::Result<PathBuf> {
    let name = format!("out_k{kmer}");
    let mut cmd = std::process::Command::new(velvetg);
    cmd.current_dir(output)
        .arg(&name)
        .arg("-cov_cutoff")
        .arg("auto")
        .arg("-exp_cov")
        .arg("auto")
        .arg("-min_contig_lgth")
        .arg("100");
    let out = std::fs::File::create(output.join(format!("velvetg-k{kmer}.out.log")))?;
    let err = std::fs::File::create(output.join(format!("velvetg-k{kmer}.err.log")))?;
    cmd.stdout(out).stderr(err);
    let status = cmd.status().context("running velvetg")?;
    anyhow::ensure!(status.success(), "velvetg failed: {status}");
    Ok(output.join(name))
}

fn cleanup_velvet_assembly_folder(output: &Path) -> anyhow::Result<()> {
    let keep = ["contigs.fa", "stats.txt"];
    for entry in std::fs::read_dir(output)? {
        let path = entry?.path();
        let name = path.file_name().unwrap().to_string_lossy().to_string();
        if !keep.contains(&name.as_str()) {
            let _ = std::fs::remove_file(&path);
        }
    }
    Ok(())
}

fn symlink(target: &Path, link: &Path) -> anyhow::Result<()> {
    let parent = link
        .parent()
        .ok_or_else(|| anyhow::anyhow!("link has no parent: {}", link.display()))?;
    let relpth = crate::assemblo_spades_cmd::pathdiff_pub(target, parent)
        .ok_or_else(|| anyhow::anyhow!("cannot compute link target for {}", link.display()))?;
    #[cfg(unix)]
    std::os::unix::fs::symlink(relpth, link)
        .with_context(|| format!("creating symlink {}", link.display()))?;
    #[cfg(not(unix))]
    anyhow::bail!("contig symlinks are only supported on Unix");
    Ok(())
}

#[allow(clippy::too_many_arguments)]
pub fn run(
    output: &Path,
    kmer: u32,
    subfolder: &str,
    clean: bool,
    config: Option<&Path>,
    dir: Option<&Path>,
) -> anyhow::Result<()> {
    let cfg = PhyluceConfig::load()?;
    let velveth = cfg.get_user_path("binaries", "velveth").map_err(|_| {
        anyhow::anyhow!("Cannot find velveth. Ensure it is installed and in your $PATH")
    })?;
    let velvetg = cfg.get_user_path("binaries", "velvetg").map_err(|_| {
        anyhow::anyhow!("Cannot find velvetg. Ensure it is installed and in your $PATH")
    })?;

    std::fs::create_dir_all(output)?;
    let contig_dir = output.join("contigs");
    std::fs::create_dir_all(&contig_dir)?;

    let input_data = get_input_data(config, dir)?;
    for (sample, sample_input_dir) in &input_data {
        crate::cli_info!("Processing {sample}");
        let sample_dir = crate::output_path::output_file(output, sample)?;
        std::fs::create_dir_all(&sample_dir)?;

        let reads = get_input_files(sample_input_dir, subfolder)?;
        let mut assembly_dir = sample_dir.clone();
        if reads.r1.is_some() && reads.r2.is_some() {
            run_velveth(&velveth, kmer, &reads, &sample_dir)?;
            assembly_dir = run_velvetg(&velvetg, kmer, &sample_dir)?;
        }
        if clean {
            cleanup_velvet_assembly_folder(&assembly_dir)?;
        }
        let contigs_file = assembly_dir.join("contigs.fa");
        if contigs_file.is_file() {
            symlink(&contigs_file, &sample_dir.join("contigs.fasta"))?;
            let link =
                crate::output_path::output_file(&contig_dir, &format!("{sample}.contigs.fasta"))?;
            symlink(&contigs_file, &link)?;
        } else {
            anyhow::bail!("Velvet did not produce contigs for sample {sample}");
        }
    }
    Ok(())
}
