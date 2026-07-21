//! CLI wiring for `phyluce align format-concatenated-phylip-for-paml`,
//! mirroring `phyluce_align_format_concatenated_phylip_for_paml`: slice a
//! relaxed-PHYLIP alignment into RAxML-style partitions and write each as
//! its own relaxed sequential PHYLIP block.
//!
//! `--config` is a RAxML partition file, e.g.
//! `DNA, p1 = 1-373, 118732-118996`; only the coordinate lists are used
//! (the `DNA, name =` prefix is ignored, matching the Python original,
//! which just splits on `=`).

use std::io::Write as _;
use std::path::Path;

use anyhow::Context;
use phyluce_align::Alignment;
use phyluce_assembly::FastMap;

type PartitionRange = (usize, usize);
type Partition = (String, Vec<PartitionRange>);

fn parse_relaxed_phylip(text: &str) -> anyhow::Result<Alignment> {
    let mut lines = text.lines();
    let header = lines
        .next()
        .ok_or_else(|| anyhow::anyhow!("empty phylip file"))?;
    let mut parts = header.split_whitespace();
    let ntax: usize = parts
        .next()
        .ok_or_else(|| anyhow::anyhow!("bad phylip header"))?
        .parse()?;
    let nchar: usize = parts
        .next()
        .ok_or_else(|| anyhow::anyhow!("bad phylip header"))?
        .parse()?;
    let mut rows = Vec::new();
    for line in lines {
        if line.trim().is_empty() {
            continue;
        }
        let mut it = line.splitn(2, char::is_whitespace);
        let id = it.next().unwrap_or("").to_string();
        let seq = it.next().unwrap_or("").trim_start().to_string();
        rows.push(phyluce_align::AlignmentRow {
            id,
            seq: seq.into_bytes(),
        });
    }
    let alignment = Alignment { rows };
    alignment.validate()?;
    anyhow::ensure!(
        alignment.ntax() == ntax,
        "PHYLIP header declares {ntax} taxa but {} were read",
        alignment.ntax()
    );
    anyhow::ensure!(
        alignment.nchar() == nchar,
        "PHYLIP header declares {nchar} characters but {} were read",
        alignment.nchar()
    );
    Ok(alignment)
}

fn slice_columns(alignment: &Alignment, start: usize, stop: usize) -> Alignment {
    let rows = alignment
        .rows
        .iter()
        .map(|r| phyluce_align::AlignmentRow {
            id: r.id.clone(),
            seq: r.seq[start..stop].to_vec(),
        })
        .collect();
    Alignment { rows }
}

fn append_columns(alignment: &mut Alignment, other: &Alignment) {
    // O(taxa) index instead of an O(taxa) `Vec::iter().find()` per row --
    // same O(taxa^2)-per-call shape as the bug fixed in
    // `phyluce-align::concat::concatenate`; see
    // `bench_append_columns_scaling_with_taxon_count`.
    let other_index: FastMap<&str, &phyluce_align::AlignmentRow> =
        other.rows.iter().map(|r| (r.id.as_str(), r)).collect();
    for row in &mut alignment.rows {
        if let Some(o) = other_index.get(row.id.as_str()) {
            row.seq.extend_from_slice(&o.seq);
        }
    }
}

fn write_relaxed_phylip(
    out: &mut impl std::io::Write,
    alignment: &Alignment,
) -> anyhow::Result<()> {
    let id_width = if alignment.rows.is_empty() {
        2
    } else {
        alignment
            .rows
            .iter()
            .map(|r| r.id.trim().len())
            .max()
            .unwrap_or(0)
            + 2
    };
    writeln!(out, "{} {}", alignment.ntax(), alignment.nchar())?;
    for row in &alignment.rows {
        writeln!(
            out,
            "{:width$}{}",
            row.id.trim(),
            std::str::from_utf8(&row.seq)?,
            width = id_width
        )?;
    }
    Ok(())
}

fn parse_partitions(config_text: &str, nchar: usize) -> anyhow::Result<Vec<Partition>> {
    let mut partitions = Vec::new();
    for line in config_text.lines() {
        let line = line.trim();
        if line.is_empty() || line.starts_with('#') {
            continue;
        }
        let (_, rhs) = line
            .split_once('=')
            .ok_or_else(|| anyhow::anyhow!("invalid partition line {line:?}"))?;
        let mut coordinates = Vec::new();
        for segment in rhs.split(',').map(str::trim) {
            let (start, stop) = segment
                .split_once('-')
                .ok_or_else(|| anyhow::anyhow!("invalid partition range {segment:?}"))?;
            let start: usize = start
                .trim()
                .parse()
                .map_err(|_| anyhow::anyhow!("invalid partition start {start:?}"))?;
            let stop: usize = stop
                .trim()
                .parse()
                .map_err(|_| anyhow::anyhow!("invalid partition stop {stop:?}"))?;
            anyhow::ensure!(
                start >= 1 && start <= stop && stop <= nchar,
                "partition range {start}-{stop} is outside 1-{}",
                nchar
            );
            coordinates.push((start, stop));
        }
        anyhow::ensure!(
            !coordinates.is_empty(),
            "partition line has no ranges: {line:?}"
        );
        partitions.push((line.to_string(), coordinates));
    }
    Ok(partitions)
}

pub fn run(phylip_alignment: &Path, config: &Path, output: &Path) -> anyhow::Result<()> {
    let phylip_text = std::fs::read_to_string(phylip_alignment)
        .with_context(|| format!("reading PHYLIP alignment {}", phylip_alignment.display()))?;
    let aln = parse_relaxed_phylip(&phylip_text)?;
    let config_text = std::fs::read_to_string(config)
        .with_context(|| format!("reading partition config {}", config.display()))?;
    let partitions = parse_partitions(&config_text, aln.nchar())?;

    let mut out = std::fs::File::create(output)
        .with_context(|| format!("creating output file {}", output.display()))?;
    for (line, coordinates) in partitions {
        let mut new_align: Option<Alignment> = None;
        for (start, stop) in coordinates {
            let slice = slice_columns(&aln, start - 1, stop);
            match &mut new_align {
                None => new_align = Some(slice),
                Some(existing) => append_columns(existing, &slice),
            }
        }
        if let Some(partition) = new_align {
            write_relaxed_phylip(&mut out, &partition)?;
            writeln!(out)?;
            crate::cli_info!("Writing partition {line}");
        }
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    // Ad hoc benchmark: `append_columns` does an `other.rows.iter().find()`
    // (linear scan) per row in `alignment`, called once per extra range in
    // a multi-range partition -- the same O(taxa^2)-per-call shape found
    // (and fixed) in `phyluce-align::concat::concatenate`. Run with:
    //   cargo +stable test --release -p phyluce-cli --lib -- --ignored --nocapture bench_append_columns
    fn synthetic_alignment(n_taxa: usize, seq_len: usize, seed: u8) -> Alignment {
        let rows = (0..n_taxa)
            .map(|i| phyluce_align::AlignmentRow {
                id: format!("taxon_{i}"),
                seq: vec![b"ACGT"[(i + seed as usize) % 4]; seq_len],
            })
            .collect();
        Alignment { rows }
    }

    #[test]
    #[ignore]
    fn bench_append_columns_scaling_with_taxon_count() {
        for n_taxa in [100usize, 200, 400, 800] {
            let mut base = synthetic_alignment(n_taxa, 100, 0);
            let extra = synthetic_alignment(n_taxa, 100, 1);
            let start = std::time::Instant::now();
            // A partition with 20 ranges calls append_columns 19 times.
            for _ in 0..19 {
                append_columns(&mut base, &extra);
            }
            let elapsed = start.elapsed();
            eprintln!(
                "[bench] append_columns x19: {n_taxa} taxa in {:?} ({:.3} ms/call)",
                elapsed,
                elapsed.as_secs_f64() * 1000.0 / 19.0
            );
        }
    }

    #[test]
    fn rejects_phylip_dimension_mismatches() {
        assert!(parse_relaxed_phylip("2 4\na AAAA\n").is_err());
        assert!(parse_relaxed_phylip("2 4\na AAAA\nb AAA\n").is_err());
    }

    #[test]
    fn rejects_invalid_partition_coordinates() {
        assert!(parse_partitions("DNA, p1 = 0-4", 4).is_err());
        assert!(parse_partitions("DNA, p1 = 3-2", 4).is_err());
        assert!(parse_partitions("DNA, p1 = 1-5", 4).is_err());
        assert!(parse_partitions("DNA, p1 = 1", 4).is_err());
        assert!(parse_partitions("DNA, p1 = 1-2, 4-4", 4).is_ok());
    }
}
