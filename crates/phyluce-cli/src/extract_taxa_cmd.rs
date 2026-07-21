//! CLI wiring for `phyluce align extract-taxa-from-alignments`, mirroring
//! `phyluce_align_extract_taxa_from_alignments`.

use std::collections::HashSet;
use std::path::Path;

use anyhow::Context;
use phyluce_align::nexus::format_nexus;
use phyluce_align::Alignment;
use phyluce_io::write_fasta_record;

use crate::informative_sites_cmd::{find_alignment_files, load_alignment};

pub fn run(
    alignments_dir: &Path,
    output_dir: &Path,
    input_format: &str,
    output_format: &str,
    exclude: &[String],
    include: &[String],
) -> anyhow::Result<()> {
    anyhow::ensure!(
        exclude.is_empty() || include.is_empty(),
        "--exclude and --include are mutually exclusive"
    );
    anyhow::ensure!(
        matches!(output_format, "fasta" | "nexus"),
        "output format '{output_format}' is not supported (only fasta/nexus)"
    );
    crate::output_path::prepare_output_dir(output_dir)?;
    let files = find_alignment_files(alignments_dir, input_format)
        .with_context(|| format!("reading alignments directory {}", alignments_dir.display()))?;

    // Mirrors `get_all_taxon_names` + `get_samples_to_run`: with neither
    // flag, keep everyone; with --exclude, keep everyone except the
    // excluded names (computed against the union of taxa across all
    // files); with --include, keep only the named taxa.
    let keep: HashSet<String> = if !include.is_empty() {
        include.iter().cloned().collect()
    } else {
        let mut all_taxa = HashSet::new();
        for file in &files {
            let alignment = load_alignment(file, input_format)
                .with_context(|| format!("loading alignment {}", file.display()))?;
            all_taxa.extend(alignment.rows.into_iter().map(|r| r.id));
        }
        if exclude.is_empty() {
            all_taxa
        } else {
            let excluded: HashSet<&String> = exclude.iter().collect();
            all_taxa
                .into_iter()
                .filter(|t| !excluded.contains(t))
                .collect()
        }
    };

    for file in &files {
        let name = file.file_name().and_then(|n| n.to_str()).unwrap_or("");
        let stem = name.split('.').next().unwrap_or(name);
        let alignment = load_alignment(file, input_format)
            .with_context(|| format!("loading alignment {}", file.display()))?;
        let filtered = Alignment {
            rows: alignment
                .rows
                .into_iter()
                .filter(|r| keep.contains(&r.id))
                .collect(),
        };
        if filtered.ntax() > 1 {
            let ext = output_format;
            let out_path = output_dir.join(format!("{stem}.{ext}"));
            if output_format == "fasta" {
                let mut out = std::fs::File::create(&out_path)
                    .with_context(|| format!("creating output file {}", out_path.display()))?;
                for row in &filtered.rows {
                    write_fasta_record(&mut out, &row.id, std::str::from_utf8(&row.seq)?)?;
                }
            } else {
                std::fs::write(&out_path, format_nexus(&filtered))
                    .with_context(|| format!("writing nexus output to {}", out_path.display()))?;
            }
        }
        print!(".");
    }
    crate::cli_info!();
    Ok(())
}
