//! CLI wiring for `phyluce align concatenate-alignments`, mirroring
//! `phyluce_align_concatenate_alignments`.

use std::path::Path;

use phyluce_align::concat::{concatenate, format_phylip, format_sets_block};
use phyluce_align::nexus::format_nexus_with_interleave;

use crate::informative_sites_cmd::{find_alignment_files, load_alignment};

pub fn run(
    alignments_dir: &Path,
    input_format: &str,
    output_dir: &Path,
    as_nexus: bool,
    as_phylip: bool,
) -> anyhow::Result<()> {
    anyhow::ensure!(
        as_nexus != as_phylip,
        "exactly one of --nexus or --phylip is required"
    );
    std::fs::create_dir_all(output_dir)?;

    let mut files = find_alignment_files(alignments_dir, input_format)?;
    files.sort();
    let mut loaded: Vec<(String, phyluce_align::Alignment)> = Vec::with_capacity(files.len());
    for file in &files {
        let name = file
            .file_name()
            .and_then(|n| n.to_str())
            .unwrap_or("")
            .to_string();
        let alignment = load_alignment(file, input_format)?;
        loaded.push((name, alignment));
    }
    // mirrors `data.sort()` on (basename, Nexus) tuples: unique basenames
    // make this equivalent to sorting by name alone.
    loaded.sort_by(|a, b| a.0.cmp(&b.0));

    let (combined, charsets) = concatenate(&loaded);

    let output_name = output_dir
        .file_name()
        .and_then(|n| n.to_str())
        .unwrap_or("concatenated");

    if as_phylip {
        let concat_path = output_dir.join(format!("{output_name}.phylip"));
        let charset_path = output_dir.join(format!("{output_name}.charsets"));
        std::fs::write(charset_path, format_sets_block(&charsets))?;
        std::fs::write(concat_path, format_phylip(&combined))?;
    } else {
        let concat_path = output_dir.join(format!("{output_name}.nexus"));
        // `Nexus.write_nexus_data` is called directly (not via
        // `Bio.AlignIO`'s `format()`), whose own `interleave` default is
        // always `False` regardless of alignment length.
        let mut text = format_nexus_with_interleave(&combined, false);
        text.push_str(&format_sets_block(&charsets));
        std::fs::write(concat_path, text)?;
    }
    Ok(())
}
