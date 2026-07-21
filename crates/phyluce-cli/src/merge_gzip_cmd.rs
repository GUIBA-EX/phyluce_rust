//! CLI wiring for `phyluce utilities merge-multiple-gzip-files`, mirroring
//! `phyluce_utilities_merge_multiple_gzip_files`.

use std::path::Path;

use phyluce_assembly::raw_reads::get_input_files;

use crate::conf::read_ini_values;

pub fn run(config: &Path, output: &Path, section: &str, trimmed: bool) -> anyhow::Result<()> {
    std::fs::create_dir_all(output)?;

    let text = std::fs::read_to_string(config)?;
    let items = read_ini_values(&text, section)?;

    if !trimmed {
        for (name, files) in &items {
            let mut sorted_files = files.clone();
            sorted_files.sort();
            let out_path = crate::output_path::output_file(output, name)?;
            let mut out = std::fs::File::create(&out_path)?;
            for infile in &sorted_files {
                let mut input = std::fs::File::open(infile)?;
                std::io::copy(&mut input, &mut out)?;
                crate::cli_info!(
                    "Copied {} to {}",
                    Path::new(infile).file_name().unwrap().to_string_lossy(),
                    name
                );
            }
        }
        return Ok(());
    }

    // --trimmed: each configured "file" is actually a directory of reads
    // already run through adapter/quality trimming, holding R1/R2/singleton
    // files. Merge each read type across those directories into a single
    // gzip file per sample, laid out the way downstream assembly commands
    // expect it: `<output>/<name>/split-adapter-quality-trimmed/`.
    for (name, paths) in &items {
        let mut sorted_paths = paths.clone();
        sorted_paths.sort();

        let mut r1_files = Vec::new();
        let mut r2_files = Vec::new();
        let mut singleton_files = Vec::new();
        for path in &sorted_paths {
            let reads = get_input_files(Path::new(path), "")?;
            if let Some(r1) = reads.r1 {
                r1_files.push(r1);
            }
            if let Some(r2) = reads.r2 {
                r2_files.push(r2);
            }
            if let Some(s) = reads.singleton {
                singleton_files.push(s);
            }
        }

        let sample_dir =
            crate::output_path::output_file(output, name)?.join("split-adapter-quality-trimmed");
        std::fs::create_dir_all(&sample_dir)?;

        for (files, read_kind) in [
            (&r1_files, "READ1"),
            (&r2_files, "READ2"),
            (&singleton_files, "READ-singleton"),
        ] {
            // Always create the output file, even with zero inputs
            // (matching the Python original, which unconditionally opens
            // it before the inner loop): some samples legitimately have no
            // singleton reads, and downstream tooling may expect all three
            // R1/R2/singleton files to exist regardless.
            let mut sorted_files = files.clone();
            sorted_files.sort();
            let new_name = format!("{name}-{read_kind}.fastq.gz");
            let mut out = std::fs::File::create(sample_dir.join(&new_name))?;
            for infile in &sorted_files {
                let mut input = std::fs::File::open(infile)?;
                std::io::copy(&mut input, &mut out)?;
                crate::cli_info!(
                    "Copied {} to {}",
                    infile.file_name().unwrap().to_string_lossy(),
                    new_name
                );
            }
        }
    }
    Ok(())
}
