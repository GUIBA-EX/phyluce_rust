//! CLI wiring for `phyluce align screen-alignments-for-problems`,
//! mirroring `phyluce_align_screen_alignments_for_problems`.

use std::path::Path;

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
) -> anyhow::Result<()> {
    crate::output_path::prepare_output_dir(output_dir)?;
    let files = find_alignment_files(alignments_dir, input_format)?;

    let mut count = 0usize;
    for file in &files {
        let name = file.file_name().and_then(|n| n.to_str()).unwrap_or("");
        let alignment = load_alignment(file, input_format)?;

        let has_n = !do_not_screen_n && alignment.rows.iter().any(|r| has_bad_base(&r.seq, b'N'));
        let has_x = !do_not_screen_x && alignment.rows.iter().any(|r| has_bad_base(&r.seq, b'X'));

        if !has_n && !has_x {
            std::fs::copy(file, output_dir.join(name))?;
            count += 1;
        } else if has_n {
            crate::cli_warn!("Removed locus {name} due to presence of N bases");
        } else {
            crate::cli_warn!("Removed locus {name} due to presence of X bases");
        }
    }
    crate::cli_info!("Copied {count} good alignments");
    Ok(())
}
