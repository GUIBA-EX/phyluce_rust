//! CLI wiring for `phyluce align convert-one-align-to-another`, mirroring
//! `phyluce_align_convert_one_align_to_another` (fasta/nexus only;
//! The phylip/phylip-relaxed/clustal/
//! emboss/stockholm output formats aren't implemented yet).

use std::collections::{BTreeMap, BTreeSet};
use std::path::{Path, PathBuf};

use anyhow::Context;
use phyluce_align::nexus::format_nexus;
use phyluce_io::write_fasta_record;

use crate::informative_sites_cmd::{find_alignment_files, load_alignment};

pub fn run(
    alignments_dir: &Path,
    output_dir: &Path,
    input_format: &str,
    output_format: &str,
    cores: usize,
    shorten_names: bool,
    name_conf: Option<&Path>,
) -> anyhow::Result<()> {
    anyhow::ensure!(cores > 0, "--cores must be greater than zero");
    anyhow::ensure!(
        output_format == "fasta" || output_format == "nexus",
        "output format '{output_format}' is not yet supported (only fasta/nexus)"
    );
    let files = find_alignment_files(alignments_dir, input_format)
        .with_context(|| format!("reading alignments directory {}", alignments_dir.display()))?;
    anyhow::ensure!(
        !files.is_empty(),
        "There are no {input_format}-formatted alignments in {}",
        alignments_dir.display()
    );
    crate::parallel::ensure_unique_output_names(files.iter().map(|file| {
        let name = file
            .file_name()
            .and_then(|name| name.to_str())
            .unwrap_or("");
        let stem = name.split('.').next().unwrap_or(name);
        format!("{stem}.{output_format}")
    }))?;

    // Validate the complete mapping before creating output files.
    let names = if shorten_names {
        Some(build_name_map(&files, input_format, name_conf)?)
    } else {
        None // Like upstream, --name-conf alone does not enable renaming.
    };
    crate::output_path::prepare_output_dir(output_dir)?;
    let count = files.len();
    crate::parallel::try_map_ordered(files, cores, |file| {
        let name = file.file_name().and_then(|n| n.to_str()).unwrap_or("");
        let stem = name.split('.').next().unwrap_or(name);
        let mut alignment = load_alignment(&file, input_format)
            .with_context(|| format!("loading alignment {}", file.display()))?;
        if let Some(names) = &names {
            for row in &mut alignment.rows {
                row.id = names
                    .get(&row.id)
                    .with_context(|| format!("no name mapping for taxon {:?}", row.id))?
                    .clone();
            }
            alignment.validate()?;
        }
        let out_path = output_dir.join(format!("{stem}.{output_format}"));
        if output_format == "fasta" {
            let mut out = std::fs::File::create(&out_path)
                .with_context(|| format!("creating output file {}", out_path.display()))?;
            for row in &alignment.rows {
                write_fasta_record(&mut out, &row.id, std::str::from_utf8(&row.seq)?)?;
            }
        } else {
            std::fs::write(&out_path, format_nexus(&alignment))
                .with_context(|| format!("writing nexus output to {}", out_path.display()))?;
        }
        Ok(())
    })?;
    for _ in 0..count {
        print!(".");
    }
    crate::cli_info!();
    if let Some(names) = names {
        for (full, short) in names {
            tracing::info!("Taxon renamed: {full} => {short}");
        }
    }
    Ok(())
}

/// Preserve the upstream separator precedence and six-character base rule.
fn shortened_name(name: &str) -> String {
    let separator = ['-', '_', ' '].into_iter().find(|c| name.contains(*c));
    if let Some(separator) = separator {
        name.split(separator)
            .take(2)
            .map(|part| {
                let mut previous_cased = false;
                let mut result = String::new();
                for c in part.chars().take(3) {
                    if previous_cased {
                        result.extend(c.to_lowercase());
                    } else {
                        result.extend(c.to_uppercase());
                    }
                    previous_cased = c.is_lowercase() || c.is_uppercase();
                }
                result
            })
            .collect()
    } else {
        name.chars().take(6).collect()
    }
}

fn build_name_map(
    files: &[PathBuf],
    input_format: &str,
    name_conf: Option<&Path>,
) -> anyhow::Result<BTreeMap<String, String>> {
    let mut taxa = BTreeSet::new();
    for file in files {
        let alignment = load_alignment(file, input_format)
            .with_context(|| format!("loading alignment {}", file.display()))?;
        taxa.extend(alignment.rows.into_iter().map(|row| row.id));
    }
    let mut names = BTreeMap::new();
    if let Some(path) = name_conf {
        let text = std::fs::read_to_string(path)
            .with_context(|| format!("reading name config {}", path.display()))?;
        let sections = crate::conf::parse_ini(&text);
        let entries = sections
            .get("taxa")
            .context("name config must contain a [taxa] section")?;
        let mut used = BTreeSet::new();
        for (full, short) in entries {
            anyhow::ensure!(!short.trim().is_empty() && !short.chars().any(char::is_whitespace),
                "invalid short name for taxon {full:?}: names must be nonempty and contain no whitespace");
            anyhow::ensure!(
                names.insert(full.clone(), short.clone()).is_none(),
                "duplicate taxon {full:?} in name config"
            );
            anyhow::ensure!(
                used.insert(short.clone()),
                "duplicate short name {short:?} in name config"
            );
        }
        for taxon in taxa {
            anyhow::ensure!(
                names.contains_key(&taxon),
                "no name mapping for taxon {taxon:?}"
            );
        }
    } else {
        let mut used = BTreeSet::new();
        for taxon in taxa {
            let base = shortened_name(&taxon);
            anyhow::ensure!(
                !base.is_empty(),
                "cannot shorten taxon {taxon:?} to an empty name"
            );
            let mut short = base.clone();
            let mut suffix = 0;
            while !used.insert(short.clone()) {
                short = format!("{base}{suffix}");
                suffix += 1;
            }
            names.insert(taxon, short);
        }
    }
    Ok(names)
}
