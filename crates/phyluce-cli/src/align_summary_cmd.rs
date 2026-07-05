//! CLI wiring for `phyluce align get-align-summary-data`, mirroring
//! `phyluce_align_get_align_summary_data`'s `--output-stats` CSV (the log
//! output isn't reproduced byte-for-byte; it isn't covered by any golden
//! fixture).

use std::path::{Path, PathBuf};

use phyluce_align::summary::compute_align_summary;

use crate::informative_sites_cmd::{find_alignment_files, load_alignment};

pub fn run(
    alignments_dir: &Path,
    input_format: &str,
    output_stats: Option<PathBuf>,
) -> anyhow::Result<()> {
    let files = find_alignment_files(alignments_dir, input_format)?;
    let mut rows = Vec::new();
    for file in &files {
        let name = file
            .file_name()
            .and_then(|n| n.to_str())
            .unwrap_or("")
            .to_string();
        let alignment = load_alignment(file, input_format)?;
        let s = compute_align_summary(&alignment);
        rows.push((name, s));
    }

    if let Some(out_path) = output_stats {
        let mut out = String::from(
            "aln,length,sites,differences,characters,gc content,gaps,a count, c count, g count, t count\n",
        );
        for (name, s) in &rows {
            out.push_str(&format!(
                "{},{},{},{},{},{},{},{},{},{},{}\n",
                name,
                s.length,
                s.sum_informative_sites,
                s.sum_differences,
                s.sum_counted_sites,
                format_gc(s.gc_content_percent()),
                s.char_count(b'-'),
                s.char_count(b'A'),
                s.char_count(b'C'),
                s.char_count(b'G'),
                s.char_count(b'T'),
            ));
        }
        std::fs::write(out_path, out)?;
    }
    Ok(())
}

/// Mirrors Python's default float `str()` formatting used by `"{}".format`
/// on a `round(x, 2)` result: whole numbers print without a trailing `.0`
/// stripped (e.g. `50.0`, not `50`), matching `round()`'s float return type.
fn format_gc(x: f64) -> String {
    if x == x.trunc() {
        format!("{x:.1}")
    } else {
        format!("{x}")
    }
}
