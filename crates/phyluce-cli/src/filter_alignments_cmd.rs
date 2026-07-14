//! CLI wiring for `phyluce align filter-alignments`, mirroring
//! `phyluce_align_filter_alignments` (copies files meeting the criteria
//! from input to output, unmodified).

use std::path::Path;

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
) -> anyhow::Result<()> {
    crate::output_path::prepare_output_dir(output_dir)?;
    let files = find_alignment_files(alignments_dir, input_format)?;

    for file in &files {
        let alignment = load_alignment(file, input_format)?;

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
            std::fs::copy(file, output_dir.join(name))?;
            crate::cli_info!("Good alignment: {}", name.to_string_lossy());
        }
    }
    Ok(())
}
