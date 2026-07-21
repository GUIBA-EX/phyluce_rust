//! CLI wiring for `phyluce align screen-alignments-for-problems`,
//! mirroring `phyluce_align_screen_alignments_for_problems`.

use std::path::Path;

use anyhow::Context;

use crate::informative_sites_cmd::{find_alignment_files, load_alignment};

fn has_bad_base(seq: &[u8], base: u8) -> bool {
    seq.iter().any(|&b| b.eq_ignore_ascii_case(&base))
}

pub fn run(
    alignments_dir: &Path,
    output_dir: &Path,
    do_not_screen_n: bool,
    do_not_screen_x: bool,
    input_format: &str,
    cores: usize,
) -> anyhow::Result<()> {
    anyhow::ensure!(cores > 0, "--cores must be greater than zero");
    crate::output_path::prepare_output_dir(output_dir)?;
    let files = find_alignment_files(alignments_dir, input_format)
        .with_context(|| format!("reading alignments directory {}", alignments_dir.display()))?;

    let results = crate::parallel::try_map_ordered(files, cores, |file| {
        let name = file.file_name().and_then(|n| n.to_str()).unwrap_or("");
        let alignment = load_alignment(&file, input_format)
            .with_context(|| format!("loading alignment {}", file.display()))?;

        let has_n = !do_not_screen_n && alignment.rows.iter().any(|r| has_bad_base(&r.seq, b'N'));
        let has_x = !do_not_screen_x && alignment.rows.iter().any(|r| has_bad_base(&r.seq, b'X'));

        if !has_n && !has_x {
            let dest = output_dir.join(name);
            std::fs::copy(&file, &dest)
                .with_context(|| format!("copying {} to {}", file.display(), dest.display()))?;
            Ok((true, None))
        } else if has_n {
            Ok((
                false,
                Some(format!("Removed locus {name} due to presence of N bases")),
            ))
        } else {
            Ok((
                false,
                Some(format!("Removed locus {name} due to presence of X bases")),
            ))
        }
    })?;
    for (_, warning) in &results {
        if let Some(warning) = warning {
            crate::cli_warn!("{warning}");
        }
    }
    let count = results.iter().filter(|(copied, _)| *copied).count();
    crate::cli_info!("Copied {count} good alignments");
    Ok(())
}
