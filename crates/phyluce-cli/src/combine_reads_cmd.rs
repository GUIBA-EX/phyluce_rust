//! CLI wiring for `phyluce utilities combine-reads`, mirroring
//! `phyluce_utilities_combine_reads`.

use std::collections::HashMap;
use std::path::Path;

use anyhow::Context;
use regex::Regex;

fn read_key(basename: &str, re: &Regex) -> Option<String> {
    let caps = re.captures(basename)?;
    if let Some(m) = caps.get(1) {
        return Some(m.as_str().to_string());
    }
    caps.get(2).map(|m| m.as_str().to_string())
}

pub fn run(config: &Path, output: &Path, subfolder: &str) -> anyhow::Result<()> {
    let text = std::fs::read_to_string(config)
        .with_context(|| format!("reading config file {}", config.display()))?;
    let directories = crate::conf::read_ini_values(&text, "samples")?;

    // mirrors `(?:.*)[_-](?:READ|Read|R)(\d)*[_-]*(singleton)*(?:.*)`
    let re = Regex::new(r"(?:.*)[_-](?:READ|Read|R)(\d)*[_-]*(singleton)*(?:.*)")?;

    for (name, locations) in &directories {
        let mut all_files: HashMap<String, Vec<std::path::PathBuf>> = HashMap::new();
        for location in locations {
            let dir = Path::new(location).join(subfolder);
            let Ok(entries) = std::fs::read_dir(&dir) else {
                continue;
            };
            for entry in entries.flatten() {
                let path = entry.path();
                let Some(fname) = path.file_name().and_then(|n| n.to_str()) else {
                    continue;
                };
                if !fname.ends_with(".fastq.gz") {
                    continue;
                }
                if let Some(key) = read_key(fname, &re) {
                    all_files.entry(key).or_default().push(path);
                }
            }
        }

        let sample_dir = crate::output_path::output_file(output, name)?;
        let output_dir = if subfolder.is_empty() {
            sample_dir
        } else {
            crate::output_path::output_file(&sample_dir, subfolder)?
        };
        std::fs::create_dir_all(&output_dir)
            .with_context(|| format!("creating output directory {}", output_dir.display()))?;
        crate::cli_info!("Processing {name}");
        for (read, files) in &all_files {
            crate::cli_info!("\tProcessing read {read}");
            let new_file_name = if read == "1" || read == "2" {
                format!("{name}-READ{read}.fastq.gz")
            } else {
                format!("{name}-READ-{read}.fastq.gz")
            };
            let new_file_path = output_dir.join(&new_file_name);
            let mut out = std::fs::File::create(&new_file_path)
                .with_context(|| format!("creating output file {}", new_file_path.display()))?;
            for (i, file) in files.iter().enumerate() {
                let mut input = std::fs::File::open(file)
                    .with_context(|| format!("opening input file {}", file.display()))?;
                std::io::copy(&mut input, &mut out)?;
                if i == 0 {
                    crate::cli_info!(
                        "\t\tCopying file {} to {}",
                        file.file_name().unwrap().to_string_lossy(),
                        new_file_name
                    );
                } else {
                    crate::cli_info!(
                        "\t\tAppending file {} to {}",
                        file.file_name().unwrap().to_string_lossy(),
                        new_file_name
                    );
                }
            }
        }
    }
    Ok(())
}
