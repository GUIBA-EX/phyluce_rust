//! Shared helper for invoking the `probebwa` binary, mirroring the
//! genome-index / hash-table / mapping workflow that `docs/tutorials/
//! tutorial-4.rst` runs by hand against `stampy.py`
//! (`stampy.py --species=... --assembly=... -G ...`, `stampy.py -g ... -H
//! ...`, `stampy.py -g ... -h ... --substitutionrate=... -M ...`).
//!
//! `probebwa` (<https://github.com/GUIBA-EX/probebwa>) is a from-scratch
//! Rust reimplementation of the published Stampy algorithm (Lunter &
//! Goodson 2011) with a deliberately stampy-compatible CLI: `build-genome`
//! / `build-hash` / `map`, with `--gapopen`/`--gapextend` even keeping
//! stampy's flag names and default values. It writes SAM/BAM directly, so
//! `easy-stampy` doesn't need a separate `samtools view` step the way the
//! tutorial's hand-run pipeline does.
//!
//! Pure argument construction, factored out so it's testable without a
//! `probebwa` binary present (same rationale as `lastz_align`).

/// `probebwa build-genome --species SPECIES --assembly ASSEMBLY -G OUTPUT FILES...`
pub fn build_genome_args(
    species: &str,
    assembly: &str,
    output_prefix: &str,
    files: &[String],
) -> Vec<String> {
    let mut args = vec![
        "build-genome".to_string(),
        "--species".to_string(),
        species.to_string(),
        "--assembly".to_string(),
        assembly.to_string(),
        "-G".to_string(),
        output_prefix.to_string(),
    ];
    args.extend(files.iter().cloned());
    args
}

/// `probebwa build-hash --genome GENOME -H OUTPUT`
pub fn build_hash_args(genome_prefix: &str, output_prefix: &str) -> Vec<String> {
    vec![
        "build-hash".to_string(),
        "--genome".to_string(),
        genome_prefix.to_string(),
        "-H".to_string(),
        output_prefix.to_string(),
    ]
}

/// `probebwa map --genome GENOME --hash HASH --substitution-rate RATE
/// --threads THREADS [--outputformat bam --output OUTPUT] -M READS...`
#[allow(clippy::too_many_arguments)]
pub fn map_args(
    genome_prefix: &str,
    hash_prefix: &str,
    substitution_rate: f64,
    threads: usize,
    reads: &[String],
    output: &str,
    bam: bool,
) -> Vec<String> {
    let mut args = vec![
        "map".to_string(),
        "--genome".to_string(),
        genome_prefix.to_string(),
        "--hash".to_string(),
        hash_prefix.to_string(),
        "--substitution-rate".to_string(),
        substitution_rate.to_string(),
        "--threads".to_string(),
        threads.to_string(),
    ];
    if bam {
        args.push("--outputformat".to_string());
        args.push("bam".to_string());
    }
    args.push("--output".to_string());
    args.push(output.to_string());
    args.push("-M".to_string());
    args.extend(reads.iter().cloned());
    args
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn builds_genome_index_args() {
        let args = build_genome_args(
            "tribolium-castaneum",
            "triCas1",
            "base/triCas1",
            &["triCas1.fasta".to_string()],
        );
        assert_eq!(
            args,
            vec![
                "build-genome",
                "--species",
                "tribolium-castaneum",
                "--assembly",
                "triCas1",
                "-G",
                "base/triCas1",
                "triCas1.fasta",
            ]
        );
    }

    #[test]
    fn builds_hash_table_args() {
        let args = build_hash_args("base/triCas1", "base/triCas1");
        assert_eq!(
            args,
            vec![
                "build-hash",
                "--genome",
                "base/triCas1",
                "-H",
                "base/triCas1"
            ]
        );
    }

    #[test]
    fn builds_paired_map_args_with_bam_output() {
        let args = map_args(
            "base/triCas1",
            "base/triCas1",
            0.05,
            8,
            &["r1.fq.gz".to_string(), "r2.fq.gz".to_string()],
            "out.bam",
            true,
        );
        assert_eq!(
            args,
            vec![
                "map",
                "--genome",
                "base/triCas1",
                "--hash",
                "base/triCas1",
                "--substitution-rate",
                "0.05",
                "--threads",
                "8",
                "--outputformat",
                "bam",
                "--output",
                "out.bam",
                "-M",
                "r1.fq.gz",
                "r2.fq.gz",
            ]
        );
    }

    #[test]
    fn builds_single_end_map_args_with_sam_output() {
        let args = map_args(
            "hg38",
            "hg38",
            0.001,
            1,
            &["probes.fq.gz".to_string()],
            "output.sam",
            false,
        );
        assert_eq!(
            args,
            vec![
                "map",
                "--genome",
                "hg38",
                "--hash",
                "hg38",
                "--substitution-rate",
                "0.001",
                "--threads",
                "1",
                "--output",
                "output.sam",
                "-M",
                "probes.fq.gz",
            ]
        );
    }
}
