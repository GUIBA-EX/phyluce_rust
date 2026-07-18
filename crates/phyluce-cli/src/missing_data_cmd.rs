//! CLI wiring for `phyluce align add-missing-data-designators`, mirroring
//! `phyluce_align_add_missing_data_designators`.

use std::collections::HashMap;
use std::path::{Path, PathBuf};

use anyhow::Context;
use phyluce_align::nexus::format_nexus;
use phyluce_align::{Alignment, AlignmentRow};

use crate::informative_sites_cmd::{find_alignment_files, load_alignment};

/// A `[section]`-per-line conf file, `allow_no_value`-style (same shape as
/// `phyluce-assembly`'s taxon-list-config parser; kept local here to avoid
/// a cross-crate dependency for a handful of lines).
fn read_section_list_config(path: &Path) -> std::io::Result<HashMap<String, Vec<String>>> {
    let text = std::fs::read_to_string(path)?;
    Ok(crate::conf::parse_ini(&text)
        .into_iter()
        .map(|(section, entries)| (section, entries.into_iter().map(|(key, _)| key).collect()))
        .collect())
}

/// Mirrors `seq.name.lstrip("_R_")`: Python's `str.lstrip` treats its
/// argument as a *character set*, not a prefix -- it strips any leading
/// run of `_` or `R` characters, not just a literal `"_R_"` prefix.
fn lstrip_r_chars(name: &str) -> &str {
    name.trim_start_matches(['_', 'R'])
}

#[allow(clippy::too_many_arguments)]
pub fn run(
    alignments_dir: &Path,
    output_dir: &Path,
    match_count_output: &Path,
    incomplete_matrix: Option<PathBuf>,
    min_taxa: usize,
    missing_character: char,
    verbatim: bool,
    input_format: &str,
    check_missing: bool,
) -> anyhow::Result<()> {
    crate::output_path::prepare_output_dir(output_dir)?;

    let config = read_section_list_config(match_count_output)?;
    let organisms: Vec<String> = config
        .get("Organisms")
        .context("no [Organisms] section in --match-count-output")?
        .iter()
        .map(|o| o.trim_end_matches('*').to_string())
        .collect();

    let missing: Option<HashMap<String, Vec<String>>> = match &incomplete_matrix {
        Some(p) => Some(read_section_list_config(p)?),
        None => None,
    };

    let files = find_alignment_files(alignments_dir, input_format)?;
    let mut dropped = Vec::new();

    for file in &files {
        let file_name = file.file_name().and_then(|n| n.to_str()).unwrap_or("");
        let locus = file_name.split('.').next().unwrap_or(file_name).to_string();
        let alignment = load_alignment(file, input_format)?;

        if alignment.ntax() < min_taxa {
            dropped.push(file_name.to_string());
            print!(".");
            continue;
        }

        let overall_length = alignment.nchar();
        let mut local_organisms = organisms.clone();
        let mut rows = Vec::with_capacity(organisms.len());
        for row in &alignment.rows {
            let stripped = lstrip_r_chars(&row.id);
            let new_name = if verbatim {
                stripped.to_string()
            } else {
                stripped.split('_').skip(1).collect::<Vec<_>>().join("_")
            };
            let pos = local_organisms
                .iter()
                .position(|o| o == &new_name)
                .with_context(|| {
                    format!("'{new_name}' (from {file_name}) is not in --match-count-output's [Organisms]")
                })?;
            local_organisms.remove(pos);
            rows.push(AlignmentRow {
                id: new_name,
                seq: row.seq.clone(),
            });
        }

        for org in &local_organisms {
            if check_missing {
                if let Some(missing_map) = &missing {
                    let has_locus = missing_map
                        .get(org)
                        .map(|loci| loci.contains(&locus))
                        .unwrap_or(false)
                        || missing_map
                            .get(&format!("{org}*"))
                            .map(|loci| loci.contains(&locus))
                            .unwrap_or(false);
                    anyhow::ensure!(
                        has_locus,
                        "locus {locus} missing from the incomplete-matrix report for organism {org}"
                    );
                }
            }
            rows.push(AlignmentRow {
                id: org.clone(),
                seq: vec![missing_character as u8; overall_length],
            });
        }

        print!(".");
        let out_path = output_dir.join(format!("{locus}.nexus"));
        std::fs::write(out_path, format_nexus(&Alignment { rows }))?;
    }
    crate::cli_info!();
    for name in &dropped {
        crate::cli_info!("Dropped {name} because of too few taxa (N < {min_taxa})");
    }
    Ok(())
}
