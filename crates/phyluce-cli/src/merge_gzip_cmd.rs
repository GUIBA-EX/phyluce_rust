//! CLI wiring for `phyluce utilities merge-multiple-gzip-files`, mirroring
//! `phyluce_utilities_merge_multiple_gzip_files`.
//!
//! Only the non-`--trimmed` path is implemented: the `--trimmed` path
//! depends on `phyluce.raw_reads`' R1/R2/singleton file-discovery logic,
//! which isn't ported yet (see docs/rust-rewrite-plan.md's phased plan).

use std::path::Path;

use crate::conf::read_ini_values;

pub fn run(config: &Path, output: &Path, section: &str, trimmed: bool) -> anyhow::Result<()> {
    anyhow::ensure!(!trimmed, "--trimmed is not yet implemented");
    std::fs::create_dir_all(output)?;

    let text = std::fs::read_to_string(config)?;
    let items = read_ini_values(&text, section)?;

    for (name, files) in &items {
        let mut sorted_files = files.clone();
        sorted_files.sort();
        let out_path = output.join(name);
        let mut out = std::fs::File::create(&out_path)?;
        for infile in &sorted_files {
            let mut input = std::fs::File::open(infile)?;
            std::io::copy(&mut input, &mut out)?;
            println!(
                "Copied {} to {}",
                Path::new(infile).file_name().unwrap().to_string_lossy(),
                name
            );
        }
    }
    Ok(())
}
