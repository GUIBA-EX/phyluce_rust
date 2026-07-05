//! CLI wiring for `phyluce ncbi prep-uce-align-files-for-ncbi`, mirroring
//! `phyluce_ncbi_prep_uce_align_files_for_ncbi` + `phyluce/ncbi.py`.
//!
//! Note: as shipped, the legacy Python command currently crashes on import
//! (`from Bio.Alphabet import IUPAC` -- `Bio.Alphabet` was removed from
//! modern Biopython), so this port has no way to golden-test against a
//! live run in this environment. Ported carefully from `phyluce/ncbi.py`'s
//! logic; treat as best-effort until validated.

use std::collections::HashMap;
use std::path::Path;

use crate::conf::parse_ini;
use crate::informative_sites_cmd::{find_alignment_files, load_alignment};

/// Bare `[section]` item list (`allow_no_value`-style), used for the
/// `[exclude taxa]`/`[exclude loci]` sections.
fn read_bare_list(sections: &HashMap<String, Vec<(String, String)>>, name: &str) -> Vec<String> {
    sections
        .get(name)
        .map(|entries| entries.iter().map(|(k, _)| k.to_lowercase()).collect())
        .unwrap_or_default()
}

/// Mirrors `configparser.ConfigParser()`'s *default* `optionxform`
/// (`str.lower`): unlike most other phyluce scripts, this one never
/// overrides it, so every key read from the conf file is lowercased
/// (values are left as-is).
fn read_kv_map(
    sections: &HashMap<String, Vec<(String, String)>>,
    name: &str,
) -> HashMap<String, String> {
    sections
        .get(name)
        .map(|entries| {
            entries
                .iter()
                .map(|(k, v)| (k.to_lowercase(), v.clone()))
                .collect()
        })
        .unwrap_or_default()
}

/// Mirrors `ncbi.get_species_name`: `(sp, species, partial, oldname)`.
fn get_species_name(id: &str, remap: &HashMap<String, String>) -> (String, String, String, String) {
    let oldname = id.to_string();
    let sp = remap.get(id).cloned().unwrap_or_else(|| id.to_string());
    // Mirrors Python's `str.capitalize()`: first char uppercased, rest
    // lowercased (not just "first letter capitalized").
    let species_raw = sp.replace('_', " ");
    let species = capitalize(&species_raw);
    let partial: String = species
        .split(' ')
        .next()
        .unwrap_or("")
        .to_lowercase()
        .chars()
        .take(3)
        .collect();
    (sp, species, partial, oldname)
}

fn capitalize(s: &str) -> String {
    let mut chars = s.chars();
    match chars.next() {
        Some(first) => first.to_uppercase().collect::<String>() + &chars.as_str().to_lowercase(),
        None => String::new(),
    }
}

/// Mirrors `ncbi.get_new_identifier`.
#[allow(clippy::too_many_arguments)]
fn get_new_identifier(
    species: &str,
    uce: &str,
    partial: &str,
    counter: usize,
    metadata: &HashMap<String, String>,
    vouchers: &HashMap<String, String>,
) -> anyhow::Result<String> {
    let title = format!("{species} ultra-conserved element locus {uce}");
    let note_template = metadata
        .get("note")
        .ok_or_else(|| anyhow::anyhow!("no 'note' in [metadata]"))?;
    // mirrors `metadata["note"].format(uce)`: replaces a bare `{}` (or
    // positional `{0}`) placeholder with the locus name.
    let note = note_template.replace("{}", uce).replace("{0}", uce);
    let moltype = metadata
        .get("moltype")
        .ok_or_else(|| anyhow::anyhow!("no 'moltype' in [metadata]"))?;
    let location = metadata
        .get("location")
        .ok_or_else(|| anyhow::anyhow!("no 'location' in [metadata]"))?;
    let voucher = vouchers
        .get(&species.to_lowercase())
        .ok_or_else(|| anyhow::anyhow!("no voucher for species '{species}'"))?;
    Ok(format!(
        "{counter}{partial} [organism={species}] [moltype={moltype}] [location={location}] [note={note}] [specimen-voucher={voucher}] {title}"
    ))
}

pub fn run(
    alignments_dir: &Path,
    conf_path: &Path,
    output: &Path,
    input_format: &str,
) -> anyhow::Result<()> {
    let conf_text = std::fs::read_to_string(conf_path)?;
    let sections = parse_ini(&conf_text);

    let remap = read_kv_map(&sections, "remap");
    let remap: HashMap<String, String> = remap
        .into_iter()
        .map(|(k, v)| (k.replace(' ', "_"), v))
        .collect();
    let metadata = read_kv_map(&sections, "metadata");
    let vouchers = read_kv_map(&sections, "vouchers");
    let taxon_excludes = read_bare_list(&sections, "exclude taxa");
    let locus_excludes = read_bare_list(&sections, "exclude loci");

    let files = find_alignment_files(alignments_dir, input_format)?;
    let mut out = std::fs::File::create(output)?;
    let mut counter = 0usize;

    for file in &files {
        let uce = file
            .file_stem()
            .and_then(|s| s.to_str())
            .unwrap_or("")
            .to_string();
        if locus_excludes.contains(&uce) {
            continue;
        }
        let alignment = load_alignment(file, input_format)?;
        for row in &alignment.rows {
            let (sp, species, partial, oldname) = get_species_name(&row.id, &remap);
            let _ = sp;
            if taxon_excludes.contains(&species) || taxon_excludes.contains(&oldname) {
                continue;
            }
            let new_id =
                get_new_identifier(&species, &uce, &partial, counter, &metadata, &vouchers)?;
            let cleaned: String = row
                .seq
                .iter()
                .filter(|&&c| c != b'-' && c != b'?')
                .map(|&c| c.to_ascii_uppercase() as char)
                .collect();
            phyluce_io::write_fasta_record(&mut out, &new_id, &cleaned)?;
            counter += 1;
        }
    }
    Ok(())
}
