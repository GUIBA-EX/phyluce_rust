//! CLI wiring for `phyluce assembly assemblo-abyss`, mirroring
//! `phyluce_assembly_assemblo_abyss`.
//!
//! Untested against live abyss-pe/abyss-se binaries in this environment
//! (not installed) -- ported carefully from source, treat as best-effort
//! until validated against a real run. `convert_abyss_contigs_to_velvet`
//! (pure FASTA post-processing, no external tool) is fully testable and
//! covered by unit tests.

use std::path::{Path, PathBuf};

use phyluce_assembly::raw_reads::get_input_files;
use phyluce_config::PhyluceConfig;
use phyluce_io::{read_fasta, write_fasta_record};

use crate::assemblo_spades_cmd::{get_input_data_pub as get_input_data, pathdiff_pub};

const IUPAC_DEGENERATE: &[char] = &[
    'B', 'D', 'H', 'K', 'M', 'S', 'R', 'W', 'V', 'Y', 'X', 'b', 'd', 'h', 'k', 'm', 's', 'r', 'w',
    'v', 'y', 'x',
];

/// Mirrors `convert_abyss_contigs_to_velvet`: rename ABySS's ` >1 length
/// cov` header fields into `NODE_{num}_length_{len}_cov_{cov}`, replace
/// any IUPAC ambiguity code with `N`, and drop contigs <= 100bp.
pub fn convert_abyss_contigs_to_velvet(contigs_file: &Path) -> anyhow::Result<PathBuf> {
    let stem = contigs_file
        .file_stem()
        .and_then(|s| s.to_str())
        .unwrap_or("");
    let out_path = contigs_file
        .parent()
        .unwrap_or(Path::new("."))
        .join(format!("{stem}-velvet.fa"));

    let records = read_fasta(contigs_file)?;
    let mut out = std::fs::File::create(&out_path)?;
    for record in &records {
        let parts: Vec<&str> = record.description.split(' ').collect();
        // mirrors `seq.description.split(" ")[:3]`: description's first
        // token is the record's own id, followed by ABySS's length/cov.
        if parts.len() < 3 {
            continue;
        }
        let (num, len, cov) = (parts[0], parts[1], parts[2]);
        let new_id = format!("NODE_{num}_length_{len}_cov_{cov}");
        let mut seq = record.sequence.clone();
        for &degen in IUPAC_DEGENERATE {
            if seq.contains(degen) {
                seq = seq.replace(degen, "N");
            }
        }
        if seq.len() > 100 {
            write_fasta_record(&mut out, &new_id, &seq)?;
        }
    }
    Ok(out_path)
}

fn run_abyss_pe(
    abyss_pe: &str,
    kmer: u32,
    reads: &phyluce_assembly::raw_reads::ReadSet,
    cores: u32,
    output: &Path,
) {
    let name = format!("out_k{kmer}");
    let mut cmd = std::process::Command::new(abyss_pe);
    cmd.current_dir(output)
        .arg(format!("k={kmer}"))
        .arg(format!("j={cores}"))
        .arg(format!("name={name}"))
        .arg(format!(
            "in={} {}",
            reads.r1.as_ref().unwrap().display(),
            reads.r2.as_ref().unwrap().display()
        ));
    if let Some(s) = &reads.singleton {
        cmd.arg(format!("se={}", s.display()));
    }
    if let (Ok(out), Ok(err)) = (
        std::fs::File::create(output.join(format!("abyss-k{kmer}.out.log"))),
        std::fs::File::create(output.join(format!("abyss-k{kmer}.err.log"))),
    ) {
        cmd.stdout(out).stderr(err);
    }
    let _ = cmd.status();
}

fn run_abyss_se(
    abyss: &str,
    kmer: u32,
    reads: &phyluce_assembly::raw_reads::ReadSet,
    output: &Path,
    abyss_se: bool,
) {
    let mut cmd = std::process::Command::new(abyss);
    cmd.current_dir(output)
        .arg("-k")
        .arg(kmer.to_string())
        .arg("-o")
        .arg(format!("out_k{kmer}-contigs.fa"))
        .arg(reads.r1.as_ref().unwrap());
    if abyss_se {
        cmd.arg(reads.r2.as_ref().unwrap());
        if let Some(s) = &reads.singleton {
            cmd.arg(s);
        }
    }
    if let (Ok(out), Ok(err)) = (
        std::fs::File::create(output.join(format!("abyss-k{kmer}.out.log"))),
        std::fs::File::create(output.join(format!("abyss-k{kmer}.err.log"))),
    ) {
        cmd.stdout(out).stderr(err);
    }
    let _ = cmd.status();
}

fn cleanup_abyss_assembly_folder(output: &Path, single_end: bool) -> anyhow::Result<()> {
    let mut keep = vec!["coverage.hist".to_string()];
    for entry in std::fs::read_dir(output)? {
        let path = entry?.path();
        let name = path.file_name().unwrap().to_string_lossy().to_string();
        if name.ends_with(".log") || name.ends_with("-stats") {
            keep.push(name);
        }
    }
    for entry in std::fs::read_dir(output)? {
        let path = entry?.path();
        if path.symlink_metadata()?.file_type().is_symlink()
            && path.extension().and_then(|e| e.to_str()) == Some("fa")
        {
            let real = std::fs::canonicalize(&path)?;
            std::fs::remove_file(&path)?;
            std::fs::rename(&real, &path)?;
            keep.push(path.file_name().unwrap().to_string_lossy().to_string());
        }
    }
    if single_end {
        for entry in std::fs::read_dir(output)? {
            let path = entry?.path();
            if path.extension().and_then(|e| e.to_str()) == Some("fa") {
                keep.push(path.file_name().unwrap().to_string_lossy().to_string());
            }
        }
    }
    for entry in std::fs::read_dir(output)? {
        let path = entry?.path();
        let name = path.file_name().unwrap().to_string_lossy().to_string();
        if !keep.contains(&name) {
            let _ = std::fs::remove_file(&path);
        }
    }
    Ok(())
}

fn find_contigs_file(output: &Path) -> anyhow::Result<PathBuf> {
    std::fs::read_dir(output)?
        .filter_map(|e| e.ok())
        .map(|e| e.path())
        .find(|p| {
            p.file_name()
                .and_then(|n| n.to_str())
                .map(|n| n.ends_with("-contigs.fa"))
                .unwrap_or(false)
        })
        .ok_or_else(|| anyhow::anyhow!("no *-contigs.fa found in {}", output.display()))
}

fn symlink(target: &Path, link: &Path) {
    #[cfg(unix)]
    let _ = std::os::unix::fs::symlink(target, link);
    #[cfg(not(unix))]
    {
        let _ = target;
        let _ = link;
    }
}

#[allow(clippy::too_many_arguments)]
pub fn run(
    output: &Path,
    kmer: u32,
    cores: u32,
    subfolder: &str,
    clean: bool,
    abyss_se: bool,
    config: Option<&Path>,
    dir: Option<&Path>,
) -> anyhow::Result<()> {
    let cfg = PhyluceConfig::load()?;
    let abyss_pe_bin = cfg.get_user_path("binaries", "abyss-pe").ok();
    let abyss_bin = cfg.get_user_path("binaries", "abyss").ok();

    std::fs::create_dir_all(output)?;
    let contig_dir = output.join("contigs");
    std::fs::create_dir_all(&contig_dir)?;

    let input_data = get_input_data(config, dir)?;
    for (sample, sample_input_dir) in &input_data {
        println!("Processing {sample}");
        let sample_dir = output.join(sample);
        std::fs::create_dir_all(&sample_dir)?;
        let reads = get_input_files(sample_input_dir, subfolder)?;

        let single_end;
        if !abyss_se && reads.r1.is_some() && reads.r2.is_some() {
            let bin = abyss_pe_bin
                .as_deref()
                .ok_or_else(|| anyhow::anyhow!("Cannot find abyss-pe"))?;
            run_abyss_pe(bin, kmer, &reads, cores, &sample_dir);
            single_end = false;
            if clean {
                cleanup_abyss_assembly_folder(&sample_dir, single_end)?;
            }
        } else {
            let bin = abyss_bin
                .as_deref()
                .ok_or_else(|| anyhow::anyhow!("Cannot find abyss"))?;
            run_abyss_se(bin, kmer, &reads, &sample_dir, abyss_se);
            single_end = true;
            if clean {
                cleanup_abyss_assembly_folder(&sample_dir, single_end)?;
            }
        }

        let contigs_file = find_contigs_file(&sample_dir)?;
        let velvet_style = convert_abyss_contigs_to_velvet(&contigs_file)?;
        symlink(&velvet_style, &sample_dir.join("contigs.fasta"));
        if let Some(relpth) = pathdiff_pub(&velvet_style, &contig_dir) {
            symlink(&relpth, &contig_dir.join(format!("{sample}.contigs.fasta")));
        }
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn converts_abyss_headers_and_masks_degenerate_bases() {
        let dir =
            std::env::temp_dir().join(format!("phyluce-abyss-convert-{}", std::process::id()));
        std::fs::create_dir_all(&dir).unwrap();
        let contigs = dir.join("out_k31-contigs.fa");
        let long_seq = "ACGT".repeat(30).replace("ACGT", "ACGY"); // > 100bp, with a degenerate base
        let content = format!(
            ">1 {} 12.3\n{}\n>2 5 1.0\nACGTN\n",
            long_seq.len(),
            long_seq
        );
        std::fs::write(&contigs, content).unwrap();

        let out = convert_abyss_contigs_to_velvet(&contigs).unwrap();
        let text = std::fs::read_to_string(out).unwrap();
        assert!(text.starts_with(">NODE_1_length_"));
        assert!(!text.contains('Y'));
        // second record is <=100bp so should be dropped
        assert!(!text.contains("NODE_2"));
    }
}
