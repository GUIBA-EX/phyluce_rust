//! CLI wiring for `phyluce utilities sample-reads-from-files`, mirroring
//! `phyluce_utilities_sample_reads_from_files`.
//!
//! Shells out to `seqkit sample` rather than the Python original's
//! `seqtk sample`: seqkit is a single static Go binary (no C toolchain/
//! htslib to build), more actively maintained, and its `sample` has the
//! same shape (a `-s` seed plus a target size). Not a byte-for-byte
//! substitute, though -- seqkit's sampling algorithm/RNG differs from
//! seqtk's, so the same seed picks a *different* read set between the two
//! tools (each is internally reproducible on repeated runs, just not
//! cross-tool). And seqkit's own `--help` warns that `-n <count>` (an
//! exact count) loads the whole FASTQ into memory on large files; this
//! samples by `-p <proportion>` instead (computed from the target count
//! and the file's own read count, via `count_fastq_reads`), which seqkit
//! streams.

use std::path::{Path, PathBuf};
use std::process::{Command, Stdio};

/// Simple xorshift RNG, consistent with the pattern used elsewhere in
/// this port (e.g. `bootstrap_count_cmd`) -- not a reproduction of
/// Python's Mersenne Twister, just a stand-in for `random.randrange`
/// (the seed is only used to keep R1/R2 in sync within one seqkit call
/// pair, not for any downstream correctness).
struct SimpleRng(u64);

impl SimpleRng {
    fn new(seed: u64) -> Self {
        Self(seed.wrapping_mul(2685821657736338717).wrapping_add(1))
    }
    fn next_u64(&mut self) -> u64 {
        self.0 ^= self.0 << 13;
        self.0 ^= self.0 >> 7;
        self.0 ^= self.0 << 17;
        self.0
    }
    fn randrange(&mut self, upper: u64) -> u64 {
        self.next_u64() % upper
    }
}

fn count_fastq_reads(path: &Path) -> anyhow::Result<usize> {
    Ok(phyluce_io::fastq::fastq_record_count(path)?)
}

#[allow(clippy::too_many_arguments)]
fn run_seqkit(
    seqkit_bin: &str,
    frac: f64,
    total_reads: i64,
    rand: u64,
    fastq: &str,
    reads: i64,
    output_index: usize,
    output_dir: &Path,
) -> anyhow::Result<()> {
    let dir_name = output_dir
        .file_name()
        .and_then(|s| s.to_str())
        .unwrap_or("");
    let read_label = if output_index == 0 { "R1" } else { "R2" };
    let out_fname = output_dir.join(format!(
        "{dir_name}_{:.0}_{total_reads}r_L001_{read_label}_001.fastq",
        frac * 100.0
    ));
    let err_fname = output_dir.join(format!(
        "{dir_name}_{:.0}_{total_reads}r_L001_{read_label}_001.seqkit-err.txt",
        frac * 100.0
    ));

    crate::cli_warn!(
        "\tfrac:{frac}, input:{fastq}, rand:{rand}, reads:{reads}, out_fname:{}",
        out_fname.display()
    );

    // seqkit's `-n <count>` loads the whole FASTQ into memory (its own
    // --help warns against this on large files); sample by proportion
    // instead, which it streams.
    let file_reads = count_fastq_reads(Path::new(fastq))?;
    anyhow::ensure!(file_reads > 0, "{fastq} has no reads to sample from");
    let proportion = (reads as f64 / file_reads as f64).clamp(0.0, 1.0);

    let out = std::fs::OpenOptions::new()
        .create(true)
        .append(true)
        .open(&out_fname)?;
    let err = std::fs::OpenOptions::new()
        .create(true)
        .append(true)
        .open(&err_fname)?;
    let seed = rand.to_string();
    let proportion_arg = proportion.to_string();
    let status = Command::new(seqkit_bin)
        .args(["sample", "-p", &proportion_arg, "-s", &seed, fastq])
        .stdout(Stdio::from(out))
        .stderr(Stdio::from(err))
        .status()?;
    anyhow::ensure!(status.success(), "seqkit failed with status {status}");
    Ok(())
}

fn sample_reads_with_seqkit(
    seqkit_bin: &str,
    frac: f64,
    total_reads: i64,
    reads: &std::collections::HashMap<String, Vec<String>>,
    to_get: &std::collections::HashMap<String, i64>,
    output_dir: &Path,
    rng: &mut SimpleRng,
) -> anyhow::Result<()> {
    crate::cli_warn!("\tGetting reads for {frac} UCE fraction");
    let mut names: Vec<&String> = to_get.keys().collect();
    names.sort();
    for name in names {
        let count = to_get[name];
        let rand = rng.randrange(1_000_000);
        if let Some(paths) = reads.get(name) {
            for (item, path) in paths.iter().enumerate() {
                run_seqkit(
                    seqkit_bin,
                    frac,
                    total_reads,
                    rand,
                    path,
                    count,
                    item,
                    output_dir,
                )?;
            }
        }
    }
    Ok(())
}

pub fn run(conf: &Path, output: &Path) -> anyhow::Result<()> {
    crate::output_path::prepare_output_dir(output)?;
    let cfg = phyluce_config::PhyluceConfig::load()?;
    let seqkit_bin = cfg.get_user_path("binaries", "seqkit")?;

    let conf_text = std::fs::read_to_string(conf)?;
    let reads = crate::conf::read_ini_values(&conf_text, "reads")?;
    let sections = crate::conf::parse_ini(&conf_text);
    let splits = sections
        .get("splits")
        .ok_or_else(|| anyhow::anyhow!("no [splits] section in --conf"))?;
    let splits_map: std::collections::HashMap<&str, &str> = splits
        .iter()
        .map(|(k, v)| (k.as_str(), v.as_str()))
        .collect();

    let total_reads: i64 = splits_map
        .get("total_reads")
        .ok_or_else(|| anyhow::anyhow!("missing splits.total_reads"))?
        .parse()?;
    let uce_frac: Vec<f64> = splits_map
        .get("uce_reads")
        .ok_or_else(|| anyhow::anyhow!("missing splits.uce_reads"))?
        .split(',')
        .map(|s| s.trim().parse::<f64>())
        .collect::<Result<_, _>>()?;
    let mtdna_frac: f64 = splits_map
        .get("mtdna_reads")
        .ok_or_else(|| anyhow::anyhow!("missing splits.mtdna_reads"))?
        .parse()?;

    let mut rng = SimpleRng::new(std::process::id() as u64 ^ 0x9E3779B97F4A7C15);

    for &frac in &uce_frac {
        let mut to_get = std::collections::HashMap::new();
        let uce = (frac * total_reads as f64) as i64;
        to_get.insert("uce".to_string(), uce);

        let mtdna_path = reads
            .get("mtdna")
            .and_then(|v| v.first())
            .ok_or_else(|| anyhow::anyhow!("missing reads.mtdna"))?;
        let mtdna_reads = count_fastq_reads(&PathBuf::from(mtdna_path))?;
        let mtdna = (mtdna_frac * mtdna_reads as f64) as i64;
        to_get.insert("mtdna".to_string(), mtdna);

        anyhow::ensure!(
            uce + mtdna <= total_reads,
            "splits.uce_reads ({frac}) and splits.mtdna_reads ({mtdna_frac}) together request {uce} + {mtdna} = {} reads, more than splits.total_reads ({total_reads})",
            uce + mtdna
        );
        let genome = total_reads - uce - mtdna;
        to_get.insert("genome".to_string(), genome);

        crate::cli_warn!(
            "Reads:{total_reads}, UCE:{uce} - {frac} on target, mtDNA:{mtdna}, genome:{genome}"
        );
        sample_reads_with_seqkit(
            &seqkit_bin,
            frac,
            total_reads,
            &reads,
            &to_get,
            output,
            &mut rng,
        )?;
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn simple_rng_stays_in_range() {
        let mut rng = SimpleRng::new(42);
        for _ in 0..100 {
            assert!(rng.randrange(1_000_000) < 1_000_000);
        }
    }
}
