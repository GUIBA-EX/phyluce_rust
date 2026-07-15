//! CLI wiring for `phyluce align get-only-loci-with-min-taxa`, mirroring
//! `phyluce_align_get_only_loci_with_min_taxa`.

use std::path::Path;

use crate::informative_sites_cmd::find_alignment_files;

pub fn run(
    alignments_dir: &Path,
    taxa: usize,
    output_dir: &Path,
    percent: f64,
    input_format: &str,
    cores: usize,
) -> anyhow::Result<()> {
    anyhow::ensure!(cores > 0, "--cores must be greater than zero");
    crate::output_path::prepare_output_dir(output_dir)?;
    let files = find_alignment_files(alignments_dir, input_format)?;
    let min_count = (percent * taxa as f64).floor() as usize;

    let total = files.len();
    let results = crate::parallel::try_map_ordered(files, cores, |file| {
        let alignment = crate::informative_sites_cmd::load_alignment(&file, input_format)?;
        if alignment.ntax() >= min_count {
            let dest = output_dir.join(file.file_name().unwrap());
            std::fs::copy(&file, dest)?;
            Ok(true)
        } else {
            Ok(false)
        }
    })?;
    let copied = results.iter().filter(|&&copied| copied).count();

    crate::cli_info!(
        "Copied {copied} alignments of {} total containing >= {percent} proportion of taxa (n = {min_count})",
        total
    );
    Ok(())
}
