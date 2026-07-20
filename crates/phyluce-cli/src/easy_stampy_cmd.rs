//! CLI wiring for `phyluce probe easy-stampy`: replaces the hand-run
//! `stampy.py` genome-prep + mapping workflow from `docs/tutorials/
//! tutorial-4.rst` with `probebwa`, a stampy-compatible Rust mapper. Runs
//! `probebwa build-genome`, `build-hash`, then `map` in sequence, so callers
//! don't have to script the three-step pipeline (or a trailing `samtools
//! view` for BAM output -- `probebwa map --outputformat=bam` writes it
//! directly).

use std::path::Path;

use phyluce_config::PhyluceConfig;
use phyluce_external::ExternalCommand;

use crate::probebwa_align::{build_genome_args, build_hash_args, map_args};

#[allow(clippy::too_many_arguments)]
pub fn run(
    species: &str,
    assembly: &str,
    genome_files: &[String],
    index_prefix: &Path,
    reads: &[String],
    substitution_rate: f64,
    threads: usize,
    output: &Path,
    bam: bool,
) -> anyhow::Result<()> {
    anyhow::ensure!(
        !reads.is_empty() && reads.len() <= 2,
        "--reads takes one file (single-end) or two (paired-end mate1 mate2), got {}",
        reads.len()
    );

    let cfg = PhyluceConfig::load()?;
    let probebwa_bin = cfg.get_user_path("binaries", "probebwa")?;
    let index_prefix = index_prefix.to_string_lossy().into_owned();
    let output = output.to_string_lossy().into_owned();

    ExternalCommand::new(&probebwa_bin)
        .args(build_genome_args(
            species,
            assembly,
            &index_prefix,
            genome_files,
        ))
        .run()?;

    ExternalCommand::new(&probebwa_bin)
        .args(build_hash_args(&index_prefix, &index_prefix))
        .run()?;

    ExternalCommand::new(&probebwa_bin)
        .args(map_args(
            &index_prefix,
            &index_prefix,
            substitution_rate,
            threads,
            reads,
            &output,
            bam,
        ))
        .run()?;

    Ok(())
}
