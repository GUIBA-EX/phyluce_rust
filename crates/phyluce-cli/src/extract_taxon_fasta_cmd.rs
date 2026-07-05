//! CLI wiring for `phyluce align extract-taxon-fasta-from-alignments`,
//! mirroring `phyluce_align_extract_taxon_fasta_from_alignments`.

use std::io::Write;
use std::path::Path;

use phyluce_io::write_fasta_record;

use crate::informative_sites_cmd::{find_alignment_files, load_alignment};

pub fn run(
    alignments_dir: &Path,
    taxon: &str,
    output: &Path,
    input_format: &str,
) -> anyhow::Result<()> {
    let files = find_alignment_files(alignments_dir, input_format)?;
    print!("Running");
    std::io::stdout().flush().ok();
    let mut outf = std::fs::File::create(output)?;
    for file in &files {
        let name = file.file_name().and_then(|n| n.to_str()).unwrap_or("");
        let locus = name.split('.').next().unwrap_or(name);
        let alignment = load_alignment(file, input_format)?;
        for row in &alignment.rows {
            if row.id != taxon {
                continue;
            }
            let seq: String = std::str::from_utf8(&row.seq)?
                .chars()
                .filter(|&c| c != '-' && c != '?')
                .map(|c| c.to_ascii_uppercase())
                .collect();
            if !seq.is_empty() {
                write_fasta_record(&mut outf, locus, &seq)?;
                print!(".");
                std::io::stdout().flush().ok();
            }
        }
    }
    println!();
    Ok(())
}
