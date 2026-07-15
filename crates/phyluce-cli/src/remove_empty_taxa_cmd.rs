//! CLI wiring for `phyluce align remove-empty-taxa`, mirroring
//! `phyluce_align_remove_empty_taxa`.

use std::path::Path;

use phyluce_align::nexus::format_nexus;
use phyluce_align::{Alignment, AlignmentRow};
use phyluce_io::write_fasta_record;

use crate::informative_sites_cmd::{find_alignment_files, load_alignment};

pub fn run(
    alignments_dir: &Path,
    output_dir: &Path,
    input_format: &str,
    output_format: &str,
    cores: usize,
) -> anyhow::Result<()> {
    anyhow::ensure!(cores > 0, "--cores must be greater than zero");
    anyhow::ensure!(
        matches!(output_format, "fasta" | "nexus"),
        "output format '{output_format}' is not supported (only fasta/nexus)"
    );
    crate::output_path::prepare_output_dir(output_dir)?;
    let files = find_alignment_files(alignments_dir, input_format)?;
    crate::parallel::ensure_unique_output_names(files.iter().map(|file| {
        let name = file
            .file_name()
            .and_then(|name| name.to_str())
            .unwrap_or("");
        let stem = name.split('.').next().unwrap_or(name);
        format!("{stem}.{output_format}")
    }))?;

    let count = files.len();
    crate::parallel::try_map_ordered(files, cores, |file| {
        let name = file.file_name().and_then(|n| n.to_str()).unwrap_or("");
        let stem = name.split('.').next().unwrap_or(name);
        let alignment = load_alignment(&file, input_format)?;

        let rows: Vec<AlignmentRow> = alignment
            .rows
            .into_iter()
            .filter(|r| {
                let all_missing = r.seq.iter().all(|&c| c == b'?');
                let all_gap = r.seq.iter().all(|&c| c == b'-');
                !(all_missing || all_gap)
            })
            .collect();
        let new_alignment = Alignment { rows };

        let ext = output_format;
        let out_path = output_dir.join(format!("{stem}.{ext}"));
        match output_format {
            "fasta" => {
                let mut out = std::fs::File::create(out_path)?;
                for row in &new_alignment.rows {
                    write_fasta_record(&mut out, &row.id, std::str::from_utf8(&row.seq)?)?;
                }
            }
            "nexus" => {
                std::fs::write(out_path, format_nexus(&new_alignment))?;
            }
            _ => unreachable!("output format was validated above"),
        }
        Ok(())
    })?;
    for _ in 0..count {
        print!(".");
    }
    crate::cli_info!();
    Ok(())
}
