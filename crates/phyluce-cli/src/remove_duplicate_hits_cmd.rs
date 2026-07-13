//! CLI wiring for `phyluce probe remove-duplicate-hits-from-probes-using-lastz`,
//! mirroring `phyluce_probe_remove_duplicate_hits_from_probes_using_lastz`.

use std::collections::{HashMap, HashSet};
use std::path::{Path, PathBuf};

use phyluce_io::lastz::read_lastz;
use phyluce_io::{read_fasta, write_fasta_record};
use regex::Regex;

fn probe_name(header: &str, regex: &Regex) -> anyhow::Result<String> {
    regex
        .captures(header)
        .and_then(|c| c.get(1))
        .map(|m| m.as_str().to_string())
        .ok_or_else(|| anyhow::anyhow!("no regex match for {:?}", header))
}

fn get_dupes(
    lastz_file: &Path,
    regex: &Regex,
    long_format: bool,
) -> anyhow::Result<HashSet<String>> {
    let matches_list = read_lastz(lastz_file, long_format)?;
    let mut matches: HashMap<String, Vec<String>> = HashMap::new();
    for m in &matches_list {
        let target = probe_name(&m.name1, regex)?;
        let query = probe_name(&m.name2, regex)?;
        matches.entry(target).or_default().push(query);
    }
    let mut dupes = HashSet::new();
    for (k, v) in &matches {
        if v.len() > 1 {
            for i in v {
                if i != k {
                    dupes.insert(k.clone());
                    dupes.insert(i.clone());
                }
            }
        } else if &v[0] != k {
            dupes.insert(k.clone());
        }
    }
    Ok(dupes.into_iter().map(|d| d.to_lowercase()).collect())
}

fn screened_filename(path: &Path) -> PathBuf {
    let parent = path.parent().unwrap_or_else(|| Path::new(""));
    let stem = path.file_stem().and_then(|s| s.to_str()).unwrap_or("");
    let ext = path
        .extension()
        .and_then(|s| s.to_str())
        .map(|e| format!(".{e}"))
        .unwrap_or_default();
    parent.join(format!("{stem}-DUPE-SCREENED{ext}"))
}

#[allow(clippy::too_many_arguments)]
pub fn run(
    fasta: &Path,
    lastz_file: &Path,
    probe_prefix: &str,
    probe_regex: &str,
    probe_bed: Option<&Path>,
    locus_bed: Option<&Path>,
    long: bool,
) -> anyhow::Result<()> {
    let regex_str = probe_regex.replace("{}", probe_prefix);
    let regex = Regex::new(&regex_str)?;
    let dupes = get_dupes(lastz_file, &regex, long)?;

    let fasta_output = screened_filename(fasta);
    let records = read_fasta(fasta)?;
    let mut out = std::fs::File::create(&fasta_output)?;
    let mut fasta_kept = 0usize;
    for record in &records {
        let locus = probe_name(&record.id, &regex)?;
        if !dupes.contains(&locus) {
            write_fasta_record(&mut out, &record.description, &record.sequence)?;
            fasta_kept += 1;
        }
    }
    crate::cli_warn!(
        "Screened {} fasta sequences.  Filtered {} duplicates. Kept {fasta_kept}.",
        records.len(),
        dupes.len()
    );

    if let Some(probe_bed) = probe_bed {
        let output = screened_filename(probe_bed);
        let text = std::fs::read_to_string(probe_bed)?;
        let mut out = std::fs::File::create(&output)?;
        use std::io::Write as _;
        let mut kept = 0usize;
        let mut count = 0usize;
        for line in text.lines() {
            if line.starts_with("track") {
                writeln!(out, "{line}")?;
                continue;
            }
            count += 1;
            let fields: Vec<&str> = line.trim().split('\t').collect();
            let locus = probe_name(fields[3], &regex)?;
            if !dupes.contains(&locus) {
                writeln!(out, "{line}")?;
                kept += 1;
            }
        }
        crate::cli_warn!(
            "Screened {count} BED probe sequences.  Filtered {} duplicates. Kept {kept}.",
            dupes.len()
        );
    }

    if let Some(locus_bed) = locus_bed {
        let output = screened_filename(locus_bed);
        let text = std::fs::read_to_string(locus_bed)?;
        let mut out = std::fs::File::create(&output)?;
        use std::io::Write as _;
        let mut kept = 0usize;
        let mut count = 0usize;
        for line in text.lines() {
            if line.starts_with("track") {
                writeln!(out, "{line}")?;
                continue;
            }
            count += 1;
            let fields: Vec<&str> = line.trim().split('\t').collect();
            let locus = fields[3];
            if !dupes.contains(locus) {
                writeln!(out, "{line}")?;
                kept += 1;
            }
        }
        crate::cli_warn!(
            "Screened {count} BED locus sequences.  Filtered {} duplicates. Kept {kept}.",
            dupes.len()
        );
    }

    Ok(())
}
