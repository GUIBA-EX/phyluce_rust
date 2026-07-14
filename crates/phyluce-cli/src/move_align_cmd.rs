//! CLI wiring for `phyluce align move-align-by-conf-file`, mirroring
//! `phyluce_align_move_align_by_conf_file`.
//!
//! The Python original doesn't override `configparser`'s default
//! `optionxform = str.lower`, so item keys from `--conf` are lowercased;
//! reproduced here via `.to_lowercase()`.

use std::collections::HashSet;
use std::path::Path;

pub fn run(
    conf: &Path,
    alignments_dir: &Path,
    output_dir: &Path,
    sections: &[String],
    opposite: bool,
    extension: &str,
) -> anyhow::Result<()> {
    crate::output_path::prepare_output_dir(output_dir)?;
    let conf_text = std::fs::read_to_string(conf)?;
    let parsed = crate::conf::parse_ini(&conf_text);

    let section_names: Vec<String> = if sections.is_empty() {
        parsed.keys().cloned().collect()
    } else {
        sections.to_vec()
    };

    let mut items: HashSet<String> = HashSet::new();
    for section in &section_names {
        if let Some(entries) = parsed.get(section) {
            for (k, _) in entries {
                items.insert(k.to_lowercase());
            }
        }
    }

    for entry in std::fs::read_dir(alignments_dir)? {
        let path = entry?.path();
        let basename = match path.file_name().and_then(|n| n.to_str()) {
            Some(n) => n,
            None => continue,
        };
        if !basename.contains(&format!(".{extension}")) {
            continue;
        }
        let matches = if !opposite {
            items.contains(basename)
        } else {
            let stem = path
                .file_stem()
                .and_then(|s| s.to_str())
                .unwrap_or(basename);
            !items.contains(stem)
        };
        if matches {
            std::fs::copy(&path, output_dir.join(basename))?;
        }
    }
    Ok(())
}
