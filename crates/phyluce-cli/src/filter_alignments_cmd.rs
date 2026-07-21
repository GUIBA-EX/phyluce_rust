//! CLI wiring for `phyluce align filter-alignments`, mirroring
//! `phyluce_align_filter_alignments` (copies files meeting the criteria
//! from input to output, unmodified).

use std::path::Path;

use anyhow::Context;

use crate::informative_sites_cmd::{find_alignment_files, load_alignment};

fn is_empty_seq(seq: &[u8]) -> bool {
    !seq.is_empty() && (seq.iter().all(|&c| c == b'-') || seq.iter().all(|&c| c == b'?'))
}

pub fn run(
    alignments_dir: &Path,
    output_dir: &Path,
    input_format: &str,
    containing: &[String],
    min_length: usize,
    min_taxa: usize,
    cores: usize,
) -> anyhow::Result<()> {
    anyhow::ensure!(cores > 0, "--cores must be greater than zero");
    crate::output_path::prepare_output_dir(output_dir)
        .with_context(|| format!("preparing output directory {}", output_dir.display()))?;
    let files = find_alignment_files(alignments_dir, input_format)
        .with_context(|| format!("reading alignments directory {}", alignments_dir.display()))?;

    let good = crate::parallel::try_map_ordered(files, cores, |file| {
        let alignment = load_alignment(&file, input_format)
            .with_context(|| format!("loading alignment {}", file.display()))?;

        let contains_ok = if containing.is_empty() {
            true
        } else {
            alignment
                .rows
                .iter()
                .any(|r| containing.contains(&r.id) && !is_empty_seq(&r.seq))
        };
        let length_ok = if min_length == 0 {
            true
        } else {
            alignment.nchar() >= min_length
        };
        let taxa_ok = if min_taxa == 0 {
            true
        } else {
            let count = alignment
                .rows
                .iter()
                .filter(|r| !is_empty_seq(&r.seq))
                .count();
            count >= min_taxa
        };

        if contains_ok && length_ok && taxa_ok {
            let name = file.file_name().unwrap();
            let dest = output_dir.join(name);
            std::fs::copy(&file, &dest)
                .with_context(|| format!("copying {} to {}", file.display(), dest.display()))?;
            Ok(Some(name.to_string_lossy().into_owned()))
        } else {
            Ok(None)
        }
    })?;
    for name in good.into_iter().flatten() {
        crate::cli_info!("Good alignment: {name}");
    }
    Ok(())
}
