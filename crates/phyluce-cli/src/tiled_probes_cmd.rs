//! CLI wiring for `phyluce probe get-tiled-probes`, mirroring
//! `phyluce_probe_get_tiled_probes`.
//!
//! As with `tiled_probe_from_multiple_inputs_cmd`, the Python original
//! breaks an odd-coordinate-count `--two-probes` tie via
//! `random.choice([1, -1])`; this port always picks `+1` instead of a
//! random pick -- a deliberate, documented divergence.

use std::io::Write as _;
use std::path::Path;

use anyhow::Context;

pub struct TiledProbesArgs {
    pub probe_prefix: String,
    pub designer: String,
    pub design: String,
    pub length: usize,
    pub density: f64,
    pub overlap_flush_left: bool,
    pub mask: Option<f64>,
    pub remove_ambiguous: bool,
    pub remove_gc: bool,
    pub start_index: usize,
    pub two_probes: bool,
}

fn middle_overlapper(seq_len: usize, length: usize, density: f64) -> Vec<(i64, i64)> {
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
    coords
}

fn validate_tiling(length: usize, density: f64) -> anyhow::Result<()> {
    anyhow::ensure!(length > 0, "--probe-length must be greater than zero");
    anyhow::ensure!(
        density.is_finite() && density > 0.0 && density <= length as f64,
        "--tiling-density must be finite, greater than zero, and no greater than --probe-length"
    );
    Ok(())
}

/// Python 3's `round()` uses banker's rounding (round-half-to-even).
fn python_round(x: f64) -> i64 {
    let floor = x.floor();
    let diff = x - floor;
    if diff < 0.5 {
        floor as i64
    } else if diff > 0.5 {
        floor as i64 + 1
    } else if (floor as i64) % 2 == 0 {
        floor as i64
    } else {
        floor as i64 + 1
    }
}

fn left_flush_overlapper(seq_len: usize, length: usize, density: f64) -> Vec<(i64, i64)> {
    let length_f = length as f64;
    let tile_overlap = length_f - length_f / density;
    let step = if tile_overlap == 0.0 {
        0
    } else {
        python_round(tile_overlap)
    };
    let stride = length as i64 - step;
    let mut coords = Vec::new();
    let mut start = 0i64;
    while start < seq_len as i64 {
        coords.push((start, start + length as i64));
        start += stride;
    }
    coords
}

struct Probe {
    id: String,
    description: String,
    seq: String,
    chromo: String,
    global_start: i64,
    global_end: i64,
    locus: String,
}

#[allow(clippy::too_many_arguments)]
pub fn run(
    input: &Path,
    output: &Path,
    probe_bed: Option<&Path>,
    locus_bed: Option<&Path>,
    args: &TiledProbesArgs,
) -> anyhow::Result<()> {
    validate_tiling(args.length, args.density)?;
    let records = phyluce_io::read_fasta(input)
        .with_context(|| format!("reading input FASTA {}", input.display()))?;

    let mut probe_set: Vec<Vec<Probe>> = Vec::new();
    for (i, record) in records.iter().enumerate() {
        let locus_count = i + args.start_index;
        let global_coords = record.description.split('|').nth(1).ok_or_else(|| {
            anyhow::anyhow!("record '{}': missing global-coords field", record.id)
        })?;
        let (global_chromo, positions) = global_coords
            .split_once(':')
            .ok_or_else(|| anyhow::anyhow!("bad global-coords field {global_coords:?}"))?;
        let global_chromo_start = positions
            .split_once('-')
            .and_then(|(s, _)| s.parse::<f64>().ok())
            .ok_or_else(|| anyhow::anyhow!("bad positions field {positions:?}"))?
            as i64;

        let seq_bytes = record.sequence.as_bytes();
        let mut coords = if args.overlap_flush_left {
            left_flush_overlapper(seq_bytes.len(), args.length, args.density)
        } else {
            middle_overlapper(seq_bytes.len(), args.length, args.density)
        };
        coords.sort();
        if args.two_probes {
            let n = coords.len();
            if n == 0 {
                continue;
            }
            coords = if n % 2 == 0 {
                vec![coords[n / 2 - 1], coords[n / 2]]
            } else {
                let pos1 = n / 2;
                let pos2 = (pos1 + 1).min(n - 1); // see module docs re: random tie-break
                let mut pair = vec![coords[pos1], coords[pos2]];
                pair.sort();
                pair
            };
        }

        let mut probes = Vec::new();
        for (k, (start, end)) in coords.into_iter().enumerate() {
            if start < 0 || end as usize > seq_bytes.len() {
                continue;
            }
            let global_probe_start = global_chromo_start + start;
            let global_probe_end = global_chromo_start + end;
            let probe_id = format!("{}{}_p{}", args.probe_prefix, locus_count, k + 1);
            let description = format!(
                " |design:{},designer:{},probes-locus:{},probes-probe:{},probes-global-chromo:{},probes-global-start:{},probes-global-end:{},probes-local-start:{},probes-local-end:{}",
                args.design, args.designer, locus_count, k + 1, global_chromo,
                global_probe_start, global_probe_end, start, end,
            );
            let slice = &seq_bytes[start as usize..end as usize];
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
                probes.push(Probe {
                    id: probe_id,
                    description,
                    seq: String::from_utf8_lossy(&upper).to_string(),
                    chromo: global_chromo.to_string(),
                    global_start: global_probe_start,
                    global_end: global_probe_end,
                    locus: locus_count.to_string(),
                });
            }
        }
        if !probes.is_empty() {
            probe_set.push(probes);
        }
    }

    let cons_count = probe_set.len();
    let probe_count: usize = probe_set.iter().map(|p| p.len()).sum();
    crate::cli_warn!("Conserved locus count = {cons_count}");
    crate::cli_warn!("Probe Count = {probe_count}");

    let mut outp = std::fs::File::create(output)
        .with_context(|| format!("creating output file {}", output.display()))?;
    let mut outpb = match probe_bed {
        Some(path) => Some(
            std::fs::File::create(path)
                .with_context(|| format!("creating probe bed output {}", path.display()))?,
        ),
        None => None,
    };
    if let Some(f) = outpb.as_mut() {
        writeln!(f, "track name=get_tiled_probes description=\"get_tiled_probes designed probes\" useScore=1 useScore=1 itemRgb=\"On\"")?;
    }
    let mut outlb = match locus_bed {
        Some(path) => Some(
            std::fs::File::create(path)
                .with_context(|| format!("creating locus bed output {}", path.display()))?,
        ),
        None => None,
    };
    if let Some(f) = outlb.as_mut() {
        writeln!(f, "track name=get_tiled_probes_loci description=\"get_tiled_probes loci\" useScore=1 useScore=1 itemRgb=\"On\"")?;
    }

    for ps in &probe_set {
        let mut lb_coords = Vec::new();
        for probe in ps {
            lb_coords.push(probe.global_start);
            lb_coords.push(probe.global_end);
            writeln!(outp, ">{}{}\n{}", probe.id, probe.description, probe.seq)?;
            if let Some(f) = outpb.as_mut() {
                writeln!(
                    f,
                    "{}\t{}\t{}\t{}\t450\t+\t0\t0\t0,0,205",
                    probe.chromo, probe.global_start, probe.global_end, probe.id
                )?;
            }
        }
        lb_coords.sort();
        let mn = *lb_coords.first().unwrap();
        let mx = *lb_coords.last().unwrap();
        anyhow::ensure!(mx > mn, "locus bed coords collapsed to a single point");
        if let Some(f) = outlb.as_mut() {
            let last = ps.last().unwrap();
            writeln!(
                f,
                "{}\t{}\t{}\t{}{}\t450\t+\t0\t0\t0,0,205",
                last.chromo, mn, mx, args.probe_prefix, last.locus
            )?;
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
