//! CLI wiring for `phyluce utilities merge-next-seq-gzip-files`, mirroring
//! `phyluce_utilities_merge_next_seq_gzip_files`.

use std::path::Path;

/// Bare `[section]` item-list parser (`allow_no_value=True`-style),
/// matching the sample-name list this command reads.
fn read_bare_section(text: &str, section: &str) -> Vec<String> {
    let mut items = Vec::new();
    let mut current: Option<String> = None;
    for raw_line in text.lines() {
        let line = raw_line.trim();
        if line.is_empty() || line.starts_with('#') {
            continue;
        }
        if line.starts_with('[') && line.ends_with(']') {
            current = Some(line[1..line.len() - 1].trim().to_string());
            continue;
        }
        if current.as_deref() == Some(section) {
            let key = line
                .split_once(':')
                .or_else(|| line.split_once('='))
                .map(|(k, _)| k.trim())
                .unwrap_or(line);
            items.push(key.to_string());
        }
    }
    items
}

fn glob_sorted(
    input_dir: &Path,
    sample: &str,
    read: &str,
) -> anyhow::Result<Vec<std::path::PathBuf>> {
    let prefix = format!("{sample}_S");
    let suffix = format!("_{read}_");
    let mut matches: Vec<std::path::PathBuf> = std::fs::read_dir(input_dir)?
        .filter_map(|e| e.ok())
        .map(|e| e.path())
        .filter(|p| {
            p.file_name()
                .and_then(|n| n.to_str())
                .map(|n| {
                    n.starts_with(&prefix)
                        && n.contains("_L")
                        && n.contains(&suffix)
                        && n.ends_with(".fastq.gz")
                })
                .unwrap_or(false)
        })
        .collect();
    matches.sort();
    Ok(matches)
}

fn merge_read(input_dir: &Path, output: &Path, sample: &str, read: &str) -> anyhow::Result<()> {
    let files = glob_sorted(input_dir, sample, read)?;
    anyhow::ensure!(
        !files.is_empty(),
        "no {read} files found for sample {sample}"
    );
    let first_name = files[0]
        .file_name()
        .unwrap()
        .to_string_lossy()
        .replace("_L001_", "_L999_");
    let out_path = output.join(&first_name);
    let mut out = std::fs::File::create(&out_path)?;
    for infile in &files {
        let mut input = std::fs::File::open(infile)?;
        std::io::copy(&mut input, &mut out)?;
        crate::cli_info!(
            "\tCopied {} to {}",
            infile.file_name().unwrap().to_string_lossy(),
            first_name
        );
    }
    Ok(())
}

pub fn run(
    input: &Path,
    config: &Path,
    output: &Path,
    section: &str,
    single_end: bool,
) -> anyhow::Result<()> {
    std::fs::create_dir_all(output)?;
    let text = std::fs::read_to_string(config)?;
    let samples = read_bare_section(&text, section);

    for sample in &samples {
        crate::cli_info!("Sample {sample}");
        merge_read(input, output, sample, "R1")?;
        if !single_end {
            merge_read(input, output, sample, "R2")?;
        }
    }
    Ok(())
}
