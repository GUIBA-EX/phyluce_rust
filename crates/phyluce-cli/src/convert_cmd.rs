//! CLI wiring for `phyluce align convert-one-align-to-another`, mirroring
//! `phyluce_align_convert_one_align_to_another` (fasta/nexus only;
//! `--shorten-names`/`--name-conf` and the phylip/phylip-relaxed/clustal/
//! emboss/stockholm output formats aren't implemented yet).

use std::path::Path;

use phyluce_align::nexus::format_nexus;
use phyluce_io::write_fasta_record;

use crate::informative_sites_cmd::{find_alignment_files, load_alignment};

pub fn run(
    alignments_dir: &Path,
    output_dir: &Path,
    input_format: &str,
    output_format: &str,
) -> anyhow::Result<()> {
    anyhow::ensure!(
        output_format == "fasta" || output_format == "nexus",
        "output format '{output_format}' is not yet supported (only fasta/nexus)"
    );
    std::fs::create_dir_all(output_dir)?;
    let files = find_alignment_files(alignments_dir, input_format)?;
    anyhow::ensure!(
        !files.is_empty(),
        "There are no {input_format}-formatted alignments in {}",
        alignments_dir.display()
    );

    for file in &files {
        let name = file.file_name().and_then(|n| n.to_str()).unwrap_or("");
        let stem = name.split('.').next().unwrap_or(name);
        let alignment = load_alignment(file, input_format)?;
        let out_path = output_dir.join(format!("{stem}.{output_format}"));
        if output_format == "fasta" {
            let mut out = std::fs::File::create(out_path)?;
            for row in &alignment.rows {
                write_fasta_record(&mut out, &row.id, std::str::from_utf8(&row.seq)?)?;
            }
        } else {
            std::fs::write(out_path, format_nexus(&alignment))?;
        }
        print!(".");
    }
    crate::cli_info!();
    Ok(())
}
