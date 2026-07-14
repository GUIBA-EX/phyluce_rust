//! CLI wiring for `phyluce align split-concat-nexus-to-loci`, mirroring
//! `phyluce_align_split_concat_nexus_to_loci`.

use std::path::Path;

use phyluce_align::charset::{parse_charsets, slice_alignment};
use phyluce_align::nexus::{format_nexus, parse_nexus};
use phyluce_io::write_fasta_record;

pub fn run(nexus_path: &Path, output_dir: &Path, output_format: &str) -> anyhow::Result<()> {
    anyhow::ensure!(
        matches!(output_format, "fasta" | "nexus"),
        "output format '{output_format}' is not supported (only fasta/nexus)"
    );
    crate::output_path::prepare_output_dir(output_dir)?;

    let text = std::fs::read_to_string(nexus_path)?;
    let mut charsets = parse_charsets(&text);
    charsets.sort_by(|a, b| a.name.cmp(&b.name));

    let alignment = parse_nexus(&text)?;

    for c in &charsets {
        let sliced = slice_alignment(&alignment, c.start, c.stop);
        let filtered = phyluce_align::Alignment {
            rows: sliced
                .rows
                .into_iter()
                .filter(|r| {
                    let all_gap = r.seq.iter().all(|&ch| ch == b'-');
                    let all_missing = r.seq.iter().all(|&ch| ch == b'?');
                    !(all_gap || all_missing)
                })
                .collect(),
        };
        let stem = c.name.replace(".nexus", "");
        let ext = output_format;
        let out_path = output_dir.join(format!("{stem}.{ext}"));
        if output_format == "fasta" {
            let mut out = std::fs::File::create(out_path)?;
            for row in &filtered.rows {
                write_fasta_record(&mut out, &row.id, std::str::from_utf8(&row.seq)?)?;
            }
        } else {
            std::fs::write(out_path, format_nexus(&filtered))?;
        }
    }
    Ok(())
}
