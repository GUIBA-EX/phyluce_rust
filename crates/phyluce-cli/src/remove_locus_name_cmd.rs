//! CLI wiring for `phyluce align remove-locus-name-from-files`, mirroring
//! `phyluce_align_remove_locus_name_from_files` (fasta/nexus output only,
//! matching the existing `convert_cmd`/`explode_alignments_cmd`
//! convention in this port).

use std::collections::HashSet;
use std::path::Path;

use phyluce_align::nexus::format_nexus;
use phyluce_align::{Alignment, AlignmentRow};
use phyluce_io::write_fasta_record;
use regex::Regex;

use crate::informative_sites_cmd::{find_alignment_files, load_alignment};

/// Mirrors Python's `str.strip("_phased")`: strips leading/trailing
/// characters that are members of the set `{_,p,h,a,s,e,d}`, not the
/// literal substring `"_phased"`.
fn strip_charset(s: &str, set: &str) -> String {
    s.trim_matches(|c| set.contains(c)).to_string()
}

pub fn run(
    alignments_dir: &Path,
    output_dir: &Path,
    taxa: Option<usize>,
    input_format: &str,
    output_format: &str,
    cores: usize,
) -> anyhow::Result<()> {
    anyhow::ensure!(cores > 0, "--cores must be greater than zero");
    anyhow::ensure!(
        output_format == "fasta" || output_format == "nexus",
        "output format '{output_format}' is not yet supported (only fasta/nexus)"
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

    print!("Running");
    let mut all_taxa: HashSet<String> = HashSet::new();
    let count = files.len();
    let taxa_sets = crate::parallel::try_map_ordered(files, cores, |file| {
        let name = file.file_name().and_then(|n| n.to_str()).unwrap_or("");
        let stem = name.split('.').next().unwrap_or(name);
        let fname = strip_charset(stem, "_phased");
        let re = Regex::new(&format!("^(_R_)*{}_*", regex::escape(&fname)))?;

        let alignment = load_alignment(&file, input_format)?;
        let mut file_taxa = HashSet::new();
        let mut rows = Vec::with_capacity(alignment.rows.len());
        for row in alignment.rows {
            let new_name = re.replacen(&row.id, 1, "").to_string();
            file_taxa.insert(new_name.clone());
            rows.push(AlignmentRow {
                id: new_name,
                seq: row.seq,
            });
        }
        if let Some(expected) = taxa {
            anyhow::ensure!(rows.len() == expected, "Taxon names are not identical");
        }
        let new_alignment = Alignment { rows };

        let out_path = output_dir.join(format!("{stem}.{output_format}"));
        if output_format == "fasta" {
            let mut out = std::fs::File::create(out_path)?;
            for row in &new_alignment.rows {
                write_fasta_record(&mut out, &row.id, std::str::from_utf8(&row.seq)?)?;
            }
        } else {
            std::fs::write(out_path, format_nexus(&new_alignment))?;
        }
        Ok(file_taxa)
    })?;
    for file_taxa in taxa_sets {
        all_taxa.extend(file_taxa);
    }
    for _ in 0..count {
        print!(".");
    }
    crate::cli_info!();
    let mut taxa_list: Vec<&String> = all_taxa.iter().collect();
    taxa_list.sort();
    crate::cli_info!(
        "Taxon names in alignments: {}",
        taxa_list
            .iter()
            .map(|s| s.as_str())
            .collect::<Vec<_>>()
            .join(",")
    );
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn strip_charset_strips_character_set_not_literal_suffix() {
        assert_eq!(strip_charset("uce-1_phased", "_phased"), "uce-1");
        assert_eq!(strip_charset("uce-1", "_phased"), "uce-1");
    }
}
