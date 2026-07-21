//! CLI wiring for `phyluce probe easy-stampy`: replaces the hand-run
//! `stampy.py` genome-prep + mapping workflow from `docs/tutorials/
//! tutorial-4.rst` with `probebwa`, a stampy-compatible Rust mapper. Runs
//! `probebwa build-genome`, `build-hash`, then `map` in sequence, so callers
//! don't have to script the three-step pipeline (or a trailing `samtools
//! view` for BAM output -- `probebwa map --outputformat=bam` writes it
//! directly).
//!
//! The whole reason the original workflow is three separate commands is so
//! `build-genome`/`build-hash` run once and `map` runs many times against
//! the same index (batch-mapping many samples to one reference). If this
//! wrapper unconditionally rebuilt the index on every call, that reuse
//! would be lost -- expensive for chromosome-scale genomes, and wasteful
//! when nothing changed. So: skip a build step whose output file already
//! exists, unless `--force-rebuild-index` says otherwise.

use std::path::Path;
use std::time::SystemTime;

use phyluce_config::PhyluceConfig;
use phyluce_external::ExternalCommand;

use crate::probebwa_align::{build_genome_args, build_hash_args, map_args};

fn mtime(path: &Path) -> Option<SystemTime> {
    std::fs::metadata(path).ok()?.modified().ok()
}

/// An index file only counts as reusable if it exists *and* isn't older
/// than any of its declared source files -- existence alone can't tell a
/// freshly rebuilt genome from a stale index left over from a previous,
/// different `--genome-files` run at the same `--index-prefix`.
fn is_fresh(index_path: &Path, sources: &[String]) -> bool {
    let Some(index_mtime) = mtime(index_path) else {
        return false;
    };
    sources
        .iter()
        .all(|s| mtime(Path::new(s)).is_some_and(|m| m <= index_mtime))
}

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
    force_rebuild_index: bool,
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

    // probebwa's own naming convention: `build-genome -G PREFIX` writes
    // `PREFIX.stidx`, `build-hash -H PREFIX` writes `PREFIX.sthash`.
    let stidx_path = format!("{index_prefix}.stidx");
    let sthash_path = format!("{index_prefix}.sthash");
    // "Reusable" means the index exists *and* isn't older than the genome
    // files it was supposedly built from -- plain existence can't tell a
    // freshly rebuilt genome from a stale index left over from a previous
    // `--genome-files` run at the same `--index-prefix` (see `is_fresh`).
    let stidx_fresh = is_fresh(Path::new(&stidx_path), genome_files);

    if force_rebuild_index || !stidx_fresh {
        ExternalCommand::new(&probebwa_bin)
            .args(build_genome_args(
                species,
                assembly,
                &index_prefix,
                genome_files,
            ))
            .run()?;
    } else {
        crate::cli_info!("Reusing existing genome index {stidx_path}");
    }

    // A rebuilt genome index invalidates any existing hash table, even if
    // the hash file itself is still on disk from a previous run -- so the
    // hash is only fresh if it's newer than the (now possibly-just-
    // rebuilt) .stidx.
    let sthash_fresh = !force_rebuild_index
        && stidx_fresh
        && is_fresh(Path::new(&sthash_path), std::slice::from_ref(&stidx_path));
    if !sthash_fresh {
        ExternalCommand::new(&probebwa_bin)
            .args(build_hash_args(&index_prefix, &index_prefix))
            .run()?;
    } else {
        crate::cli_info!("Reusing existing hash table {sthash_path}");
    }

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

#[cfg(test)]
mod tests {
    use super::*;
    use std::time::Duration;

    fn touch(path: &Path) {
        std::fs::write(path, b"x").unwrap();
    }

    fn age(path: &Path, ago: Duration) {
        let t = SystemTime::now() - ago;
        let f = std::fs::File::open(path).unwrap();
        f.set_modified(t).unwrap();
    }

    #[test]
    fn index_missing_is_never_fresh() {
        let dir =
            std::env::temp_dir().join(format!("phyluce-stampy-fresh-{}-a", std::process::id()));
        std::fs::create_dir_all(&dir).unwrap();
        let genome = dir.join("genome.fasta");
        touch(&genome);
        let index = dir.join("idx.stidx");
        assert!(!is_fresh(&index, &[genome.to_string_lossy().into_owned()]));
    }

    #[test]
    fn index_older_than_source_is_stale() {
        let dir =
            std::env::temp_dir().join(format!("phyluce-stampy-fresh-{}-b", std::process::id()));
        std::fs::create_dir_all(&dir).unwrap();
        let index = dir.join("idx.stidx");
        touch(&index);
        age(&index, Duration::from_secs(3600));
        let genome = dir.join("genome.fasta");
        touch(&genome);
        assert!(!is_fresh(&index, &[genome.to_string_lossy().into_owned()]));
    }

    #[test]
    fn index_newer_than_source_is_fresh() {
        let dir =
            std::env::temp_dir().join(format!("phyluce-stampy-fresh-{}-c", std::process::id()));
        std::fs::create_dir_all(&dir).unwrap();
        let genome = dir.join("genome.fasta");
        touch(&genome);
        age(&genome, Duration::from_secs(3600));
        let index = dir.join("idx.stidx");
        touch(&index);
        assert!(is_fresh(&index, &[genome.to_string_lossy().into_owned()]));
    }
}
