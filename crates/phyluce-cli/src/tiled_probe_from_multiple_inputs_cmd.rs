//! CLI wiring for `phyluce probe get-tiled-probe-from-multiple-inputs`,
//! mirroring `phyluce_probe_get_tiled_probe_from_multiple_inputs`.
//!
//! The Python original breaks a tie (odd number of tiling coordinates, with
//! `--two-probes`) via `random.choice([1, -1])`; this port always picks `+1`
//! instead of a random pick -- a deliberate, documented divergence.

use std::collections::HashMap;
use std::io::Write as _;
use std::path::{Path, PathBuf};

use anyhow::Context;
use phyluce_io::FastaRecord;

pub struct TilingArgs {
    // Accepted for CLI parity with the Python original, which never
    // actually reads `args.probe_prefix` either.
    #[allow(dead_code)]
    pub probe_prefix: String,
    pub designer: String,
    pub design: String,
    pub length: usize,
    pub density: f64,
    pub mask: Option<f64>,
    pub remove_ambiguous: bool,
    pub remove_gc: bool,
    pub start_index: usize,
    pub two_probes: bool,
}

fn get_fasta_files(dir: &Path) -> anyhow::Result<Vec<PathBuf>> {
    let mut out = Vec::new();
    for entry in std::fs::read_dir(dir)
        .with_context(|| format!("reading fastas directory {}", dir.display()))?
    {
        let path = entry?.path();
        if let Some(ext) = path.extension().and_then(|e| e.to_str()) {
            if matches!(ext, "fasta" | "fsa" | "aln" | "fa") {
                out.push(path);
            }
        }
    }
    out.sort();
    Ok(out)
}

fn middle_overlapper(
    seq_len: usize,
    length: usize,
    density: f64,
) -> anyhow::Result<Vec<(i64, i64)>> {
    let seq_len = seq_len as f64;
    let length_f = length as f64;
    let tile_overlap = length_f - length_f / density;
    let tile_non_overlap = length_f - tile_overlap;
    let mut coords = Vec::new();
    let middle = seq_len / 2.0;
    let halfsies = tile_overlap / 2.0;
    let mut r_prb_strt = middle - halfsies;
    let mut l_prb_strt = middle + halfsies;
    while r_prb_strt + length_f <= seq_len {
        let end = r_prb_strt + length_f;
        coords.push((r_prb_strt as i64, end as i64));
        r_prb_strt += tile_non_overlap;
    }
    while l_prb_strt - length_f >= 0.0 {
        let start = l_prb_strt - length_f;
        coords.push((start as i64, l_prb_strt as i64));
        l_prb_strt -= tile_non_overlap;
    }
    anyhow::ensure!(
        !coords.is_empty(),
        "Ensure your tiling density is sensible."
    );
    Ok(coords)
}

fn validate_tiling(length: usize, density: f64) -> anyhow::Result<()> {
    anyhow::ensure!(length > 0, "--probe-length must be greater than zero");
    anyhow::ensure!(
        density.is_finite() && density > 0.0 && density <= length as f64,
        "--tiling-density must be finite, greater than zero, and no greater than --probe-length"
    );
    Ok(())
}

struct LocusMeta {
    chromo: String,
    start: i64,
    source: String,
    seq: Vec<u8>,
}

fn parse_locus_meta(record: &FastaRecord, taxon: &str) -> anyhow::Result<LocusMeta> {
    let fields: Vec<&str> = record.description.split('|').collect();
    anyhow::ensure!(
        fields.len() >= 3,
        "malformed locus header: {:?}",
        record.description
    );
    let chromo = fields[1]
        .split_once(':')
        .map(|(_, v)| v.to_string())
        .ok_or_else(|| anyhow::anyhow!("missing chromo field"))?;
    let coords = fields[2]
        .split_once(':')
        .map(|(_, v)| v.to_string())
        .ok_or_else(|| anyhow::anyhow!("missing coords field"))?;
    let start = coords
        .split_once('-')
        .and_then(|(s, _)| s.parse::<f64>().ok())
        .ok_or_else(|| anyhow::anyhow!("bad coords field {coords:?}"))? as i64;
    Ok(LocusMeta {
        chromo,
        start,
        source: taxon.to_string(),
        seq: record.sequence.clone().into_bytes(),
    })
}

fn design_probes(
    locus_name: &str,
    loci: &[LocusMeta],
    args: &TilingArgs,
) -> anyhow::Result<Vec<(String, String)>> {
    let mut k = args.start_index;
    let mut probes = Vec::new();
    for locus in loci {
        let mut coords = middle_overlapper(locus.seq.len(), args.length, args.density)?;
        coords.sort();
        if args.two_probes {
            let n = coords.len();
            coords = if n % 2 == 0 {
                let pos1 = n / 2 - 1;
                let pos2 = n / 2;
                vec![coords[pos1], coords[pos2]]
            } else {
                let pos1 = n / 2;
                let pos2 = (pos1 as i64 + 1) as usize; // see module docs re: random tie-break
                let mut pair = vec![coords[pos1], coords[pos2.min(n - 1)]];
                pair.sort();
                pair
            };
        }
        for (start, end) in coords {
            let global_start = locus.start + start;
            let global_end = locus.start + end;
            let probe_id = format!("{locus_name}_p{k}");
            let description = format!(
                " |design:{},designer:{},probes-locus:{},probes-probe:{},probes-source:{},probes-global-chromo:{},probes-global-start:{},probes-global-end:{},probes-local-start:{},probes-local-end:{}",
                args.design, args.designer, locus_name, k, locus.source, locus.chromo,
                global_start, global_end, start, end,
            );
            let slice: Vec<u8> = locus.seq[start as usize..end as usize].to_vec();
            let upper: Vec<u8> = slice.iter().map(|b| b.to_ascii_uppercase()).collect();
            let masked =
                slice.iter().filter(|b| b.is_ascii_lowercase()).count() as f64 / slice.len() as f64;
            let gc = upper.iter().filter(|&&b| b == b'C' || b == b'G').count() as f64
                / upper.len() as f64;
            let masked_out = args.mask.map(|m| masked >= m).unwrap_or(false);
            let ambiguous_out = args.remove_ambiguous && upper.contains(&b'N');
            let gc_out = args.remove_gc && !(0.3..=0.7).contains(&gc);
            let keep = !masked_out && !ambiguous_out && !gc_out && upper.len() >= args.length;
            if keep {
                probes.push((
                    format!("{}{}", probe_id, description),
                    String::from_utf8_lossy(&upper).to_string(),
                ));
            }
            k += 1;
        }
    }
    Ok(probes)
}

#[allow(clippy::too_many_arguments)]
pub fn run(
    fastas_dir: &Path,
    multi_fasta_output: &Path,
    output: &Path,
    args: &TilingArgs,
) -> anyhow::Result<()> {
    validate_tiling(args.length, args.density)?;
    let conf_text = std::fs::read_to_string(multi_fasta_output).with_context(|| {
        format!(
            "reading multi-fasta-output file {}",
            multi_fasta_output.display()
        )
    })?;
    let sections = crate::conf::parse_ini(&conf_text);
    let hits = sections
        .get("hits")
        .ok_or_else(|| anyhow::anyhow!("no [hits] section in --multi-fasta-output"))?;
    let mut d: HashMap<String, Vec<LocusMeta>> = hits
        .iter()
        .map(|(name, _)| (name.clone(), Vec::new()))
        .collect();

    for file in get_fasta_files(fastas_dir)? {
        let taxon_name = file
            .file_stem()
            .and_then(|s| s.to_str())
            .unwrap_or("")
            .to_string();
        let records = phyluce_io::read_fasta(&file)
            .with_context(|| format!("reading fasta {}", file.display()))?;
        for record in &records {
            let fields: Vec<&str> = record.id.split('|').collect();
            if fields.len() < 4 {
                continue;
            }
            let locus_name = fields[3]
                .split_once(':')
                .map(|(_, v)| v.to_string())
                .unwrap_or_default();
            if let Some(bucket) = d.get_mut(&locus_name) {
                bucket.push(parse_locus_meta(record, &taxon_name)?);
            }
        }
    }

    let mut probe_set: Vec<Vec<(String, String)>> = Vec::new();
    let mut names: Vec<&String> = d.keys().collect();
    names.sort();
    for locus_name in names {
        probe_set.push(design_probes(locus_name, &d[locus_name], args)?);
    }

    let cons_count = probe_set.len();
    let probe_count: usize = probe_set.iter().map(|p| p.len()).sum();
    crate::cli_warn!("Conserved locus count = {cons_count}");
    crate::cli_warn!("Probe Count = {probe_count}");

    let mut out = std::fs::File::create(output)
        .with_context(|| format!("creating output file {}", output.display()))?;
    for ps in &probe_set {
        for (header, seq) in ps {
            writeln!(out, ">{header}\n{seq}")?;
        }
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn rejects_non_progressing_tiling_parameters() {
        assert!(validate_tiling(0, 2.0).is_err());
        assert!(validate_tiling(120, 0.0).is_err());
        assert!(validate_tiling(120, f64::INFINITY).is_err());
        assert!(validate_tiling(120, 121.0).is_err());
        assert!(validate_tiling(120, 2.0).is_ok());
    }
}
