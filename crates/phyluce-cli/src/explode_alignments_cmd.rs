//! CLI wiring for `phyluce align explode-alignments`, mirroring
//! `phyluce_align_explode_alignments`: convert a directory of alignments
//! into per-locus or per-taxon FASTA files, with gaps/`?` stripped.
//!
//! Note: the Python original (`bin/align/phyluce_align_explode_alignments`
//! lines 112-118, 133-136) only builds its `names` rename map when `--conf`
//! is passed, yet unconditionally does `names[name]`/`except: shortname =
//! name` further down -- if `--conf` is omitted the `names` lookup would
//! raise `NameError` rather than the intended "fall back to raw taxon
//! name". We treat "no --conf" as "no renaming" throughout, which matches
//! the evident intent of the `except` fallback.

use std::collections::HashMap;
use std::io::Write;
use std::path::{Path, PathBuf};

use crate::informative_sites_cmd::{find_alignment_files, load_alignment};

#[allow(clippy::too_many_arguments)]
pub fn run(
    alignments_dir: &Path,
    output_dir: &Path,
    input_format: &str,
    conf: Option<PathBuf>,
    section: Option<String>,
    exclude: &[String],
    by_taxon: bool,
    include_locus: bool,
) -> anyhow::Result<()> {
    crate::output_path::prepare_output_dir(output_dir)?;
    let files = find_alignment_files(alignments_dir, input_format)?;

    let mut names: HashMap<String, String> = HashMap::new();
    if let Some(conf_path) = conf {
        let section =
            section.ok_or_else(|| anyhow::anyhow!("--section is required with --conf"))?;
        let text = std::fs::read_to_string(&conf_path)?;
        let ini = phyluce_config::Ini::parse(&text);
        if let Some(entries) = ini.entries(&section) {
            for (k, v) in entries {
                names.insert(k.replace(' ', "_"), v.clone());
            }
        }
        crate::cli_info!("Original taxon count =  {}", names.len());
        for taxon in exclude {
            names.remove(taxon);
        }
    }

    let resolve =
        |name: &str| -> String { names.get(name).cloned().unwrap_or_else(|| name.to_string()) };

    if by_taxon {
        let mut handles: HashMap<String, std::fs::File> = HashMap::new();
        for file in &files {
            print!(".");
            std::io::stdout().flush().ok();
            let name = file.file_name().and_then(|n| n.to_str()).unwrap_or("");
            let locus = name.split('.').next().unwrap_or(name).to_string();
            let alignment = load_alignment(file, input_format)?;
            for row in &alignment.rows {
                let taxon_name = row
                    .id
                    .replace(&locus, "")
                    .replace("_R_", "")
                    .trim_start_matches('_')
                    .to_string();
                if exclude.iter().any(|e| e == &taxon_name) {
                    continue;
                }
                let shortname = resolve(&taxon_name);
                let seq: String = std::str::from_utf8(&row.seq)?
                    .chars()
                    .filter(|&c| c != '-' && c != '?')
                    .collect();
                if seq.is_empty() {
                    continue;
                }
                let handle = if let Some(h) = handles.get_mut(&shortname) {
                    h
                } else {
                    let path =
                        crate::output_path::output_file(output_dir, &format!("{shortname}.fasta"))?;
                    handles.insert(shortname.clone(), std::fs::File::create(path)?);
                    handles.get_mut(&shortname).unwrap()
                };
                if include_locus {
                    writeln!(handle, ">{locus}_{0} |{locus}", row.id)?;
                } else {
                    writeln!(handle, ">{}", row.id)?;
                }
                writeln!(handle, "{seq}")?;
            }
        }
    } else {
        let mut taxon_counts: Vec<usize> = Vec::new();
        for file in &files {
            print!(".");
            std::io::stdout().flush().ok();
            let name = file.file_name().and_then(|n| n.to_str()).unwrap_or("");
            let locus = name.split('.').next().unwrap_or(name).to_string();
            let alignment = load_alignment(file, input_format)?;
            let out_path = output_dir.join(format!("{locus}.fasta"));
            let mut outp = std::fs::File::create(out_path)?;
            let mut count = 0usize;
            for row in &alignment.rows {
                let taxon_name = row
                    .id
                    .replace(&locus, "")
                    .replace("_R_", "")
                    .trim_start_matches('_')
                    .to_string();
                if exclude.iter().any(|e| e == &taxon_name) {
                    continue;
                }
                let shortname = resolve(&taxon_name);
                let seq: String = std::str::from_utf8(&row.seq)?
                    .chars()
                    .filter(|&c| c != '-' && c != '?')
                    .collect();
                if seq.is_empty() {
                    crate::cli_info!("{locus}");
                    continue;
                }
                writeln!(outp, ">{shortname}")?;
                writeln!(outp, "{seq}")?;
                count += 1;
            }
            taxon_counts.push(count);
        }
        crate::cli_info!();
        let unique: std::collections::HashSet<usize> = taxon_counts.into_iter().collect();
        crate::cli_info!("Final taxon count =  {unique:?}");
    }
    Ok(())
}
