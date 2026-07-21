//! CLI wiring for `phyluce probe get-probe-bed-from-lastz-files` and
//! `phyluce probe get-locus-bed-from-lastz-files`, mirroring
//! `phyluce_probe_get_probe_bed_from_lastz_files` /
//! `phyluce_probe_get_locus_bed_from_lastz_files`.

use std::collections::{HashMap, HashSet};
use std::io::Write as _;
use std::path::Path;

use anyhow::Context;
use phyluce_io::lastz::read_lastz;
use regex::Regex;

/// Mirrors `re.search("_v_([A-Za-z0-9]+).lastz", basename).groups()[0]`.
fn outname_from_filename(name: &str) -> Option<String> {
    let re = Regex::new(r"_v_([A-Za-z0-9]+)\.lastz").ok()?;
    re.captures(name)
        .and_then(|c| c.get(1))
        .map(|m| m.as_str().to_string())
}

fn find_lastz_files(dir: &Path) -> std::io::Result<Vec<std::path::PathBuf>> {
    let mut files: Vec<_> = std::fs::read_dir(dir)?
        .filter_map(|e| e.ok())
        .map(|e| e.path())
        .filter(|p| {
            p.file_name()
                .and_then(|n| n.to_str())
                .map(|n| n.contains("lastz"))
                .unwrap_or(false)
        })
        .collect();
    files.sort();
    Ok(files)
}

fn write_bed_line(
    out: &mut std::fs::File,
    chromo: &str,
    start: i64,
    end: i64,
    name: &str,
) -> std::io::Result<()> {
    writeln!(
        out,
        "{chromo}\t{start}\t{end}\t{name}\t1000\t+\t{start}\t{end}\t100,149,237"
    )
}

/// Mirrors `phyluce_probe_get_probe_bed_from_lastz_files`.
pub fn run_probe_bed(alignments: &Path, output: &Path) -> anyhow::Result<()> {
    std::fs::create_dir_all(output)
        .with_context(|| format!("creating output directory {}", output.display()))?;
    for file in find_lastz_files(alignments)
        .with_context(|| format!("reading lastz directory {}", alignments.display()))?
    {
        let name = file.file_name().and_then(|n| n.to_str()).unwrap_or("");
        let Some(outname) = outname_from_filename(name) else {
            continue;
        };
        crate::cli_info!("Working on {outname}");

        let matches = read_lastz(&file, true)
            .with_context(|| format!("reading lastz file {}", file.display()))?;
        let mut probes: HashMap<String, Vec<(String, i64, i64)>> = HashMap::new();
        for m in &matches {
            let probe = m.name2.split('|').next().unwrap_or("").trim().to_string();
            probes
                .entry(probe)
                .or_default()
                .push((m.name1.clone(), m.zstart1, m.end1));
        }

        let out_path = output.join(format!("{outname}.probe.bed"));
        let mut out = std::fs::File::create(&out_path)
            .with_context(|| format!("creating output file {}", out_path.display()))?;
        writeln!(
            out,
            "track name=\"uce-v-{outname}\" description=\"UCE probe matches to {outname}\" visibility=2 itemRgb=\"On\""
        )?;
        let mut written = HashSet::new();
        let mut probe_names: Vec<&String> = probes.keys().collect();
        probe_names.sort();
        for probe in probe_names {
            for (chromo, start, end) in &probes[probe] {
                if written.contains(probe) {
                    crate::cli_warn!("{probe} may have >1 hit");
                } else {
                    written.insert(probe.clone());
                }
                write_bed_line(&mut out, chromo, *start, *end, probe)?;
            }
        }
    }
    Ok(())
}

/// Mirrors `phyluce_probe_get_locus_bed_from_lastz_files`.
pub fn run_locus_bed(alignments: &Path, output: &Path, regex_str: &str) -> anyhow::Result<()> {
    std::fs::create_dir_all(output)
        .with_context(|| format!("creating output directory {}", output.display()))?;
    let regex = Regex::new(regex_str)?;

    for file in find_lastz_files(alignments)
        .with_context(|| format!("reading lastz directory {}", alignments.display()))?
    {
        let name = file.file_name().and_then(|n| n.to_str()).unwrap_or("");
        let Some(outname) = outname_from_filename(name) else {
            continue;
        };
        crate::cli_info!("Working on {outname}");

        let matches = read_lastz(&file, true)
            .with_context(|| format!("reading lastz file {}", file.display()))?;
        // locus -> chromo -> positions
        let mut loci: HashMap<String, HashMap<String, Vec<i64>>> = HashMap::new();
        for m in &matches {
            let probe_field = m.name2.split('|').next().unwrap_or("").trim();
            let locus = regex
                .captures(probe_field)
                .and_then(|c| c.get(1))
                .map(|g| g.as_str().to_string())
                .ok_or_else(|| anyhow::anyhow!("no regex match for probe field {probe_field:?}"))?;
            let entry = loci
                .entry(locus)
                .or_default()
                .entry(m.name1.clone())
                .or_default();
            entry.push(m.zstart1);
            entry.push(m.end1);
        }

        let out_path = output.join(format!("{outname}.bed"));
        let mut out = std::fs::File::create(&out_path)
            .with_context(|| format!("creating output file {}", out_path.display()))?;
        writeln!(
            out,
            "track name=\"uce-v-{outname}\" description=\"UCE locus matches to {outname}\" visibility=2 itemRgb=\"On\""
        )?;
        let mut written = HashSet::new();
        let mut locus_names: Vec<&String> = loci.keys().collect();
        locus_names.sort();
        for locus in locus_names {
            let matches_for_locus = &loci[locus];
            for (chromo, locs) in matches_for_locus {
                if written.contains(locus) {
                    crate::cli_warn!("{locus} may have >1 hit");
                } else {
                    written.insert(locus.clone());
                }
                let mn = *locs.iter().min().unwrap();
                let mx = *locs.iter().max().unwrap();
                if locs.len() != 2 && mx - mn > 1000 {
                    crate::cli_warn!(
                        "Region ({} bp) is large for {locus} at {chromo}:{mn}-{mx}",
                        mx - mn
                    );
                }
                write_bed_line(&mut out, chromo, mn, mx, locus)?;
            }
        }
    }
    Ok(())
}
