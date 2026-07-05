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

use phyluce_align::Alignment;

fn parse_relaxed_phylip(text: &str) -> anyhow::Result<Alignment> {
    let mut lines = text.lines();
    let header = lines
        .next()
        .ok_or_else(|| anyhow::anyhow!("empty phylip file"))?;
    let mut parts = header.split_whitespace();
    let _ntax: usize = parts
        .next()
        .ok_or_else(|| anyhow::anyhow!("bad phylip header"))?
        .parse()?;
    let _nchar: usize = parts
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
    Ok(Alignment { rows })
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
    for row in &mut alignment.rows {
        if let Some(o) = other.rows.iter().find(|r| r.id == row.id) {
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

pub fn run(phylip_alignment: &Path, config: &Path, output: &Path) -> anyhow::Result<()> {
    let aln = parse_relaxed_phylip(&std::fs::read_to_string(phylip_alignment)?)?;
    let config_text = std::fs::read_to_string(config)?;

    let mut out = std::fs::File::create(output)?;
    for line in config_text.lines() {
        let line = line.trim();
        if line.is_empty() {
            continue;
        }
        let Some((_, rhs)) = line.split_once('=') else {
            continue;
        };
        let segments: Vec<&str> = rhs.split(',').map(|s| s.trim()).collect();
        let coordinates: Vec<(usize, usize)> = segments
            .iter()
            .map(|seg| {
                let mut it = seg.split('-');
                let a: usize = it.next().unwrap_or("1").parse().unwrap_or(1);
                let b: usize = it.next().unwrap_or("1").parse().unwrap_or(1);
                (a, b)
            })
            .collect();

        let mut new_align: Option<Alignment> = None;
        for (a, b) in &coordinates {
            let slice = slice_columns(&aln, a - 1, *b);
            match &mut new_align {
                None => new_align = Some(slice),
                Some(existing) => append_columns(existing, &slice),
            }
        }
        if let Some(partition) = new_align {
            write_relaxed_phylip(&mut out, &partition)?;
            writeln!(out)?;
            println!("Writing partition {line}");
        }
    }
    Ok(())
}
