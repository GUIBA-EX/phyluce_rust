//! CLI wiring for `phyluce utilities unmix-fasta-reads`, mirroring
//! `phyluce_utilities_unmix_fasta_reads`.

use std::io::Write as _;
use std::path::Path;

use anyhow::Context;
use phyluce_io::FastaRecord;

fn read_fasta_records(path: &Path) -> anyhow::Result<Vec<FastaRecord>> {
    phyluce_io::read_fasta(path).with_context(|| format!("reading fasta {}", path.display()))
}

#[allow(clippy::too_many_arguments)]
pub fn run(
    mixed_reads: &Path,
    singleton_reads: Option<&Path>,
    out_r1: &Path,
    out_r2: &Path,
    out_r_singleton: &Path,
    new_style: bool,
) -> anyhow::Result<()> {
    let mut records = read_fasta_records(mixed_reads)?;

    if new_style {
        for r in &mut records {
            let dsplit: Vec<&str> = r.description.split(' ').collect();
            if let Some(second) = dsplit.get(1) {
                let readnum = second.split(':').next().unwrap_or("");
                let new_id = format!("{}/{}", r.id, readnum);
                let rest = dsplit[1..].join(" ");
                r.description = format!("{new_id} {rest}");
                r.id = new_id;
            }
        }
    }

    // sort by id, matching `sorted(seq_dict.keys())`
    records.sort_by(|a, b| a.id.cmp(&b.id));

    let mut out_r1_f = std::fs::File::create(out_r1)
        .with_context(|| format!("creating R1 output file {}", out_r1.display()))?;
    let mut out_r2_f = std::fs::File::create(out_r2)
        .with_context(|| format!("creating R2 output file {}", out_r2.display()))?;
    let mut out_rs_f = std::fs::File::create(out_r_singleton).with_context(|| {
        format!(
            "creating singleton output file {}",
            out_r_singleton.display()
        )
    })?;

    if let Some(singleton) = singleton_reads {
        let contents = std::fs::read_to_string(singleton)
            .with_context(|| format!("reading singleton reads {}", singleton.display()))?;
        out_rs_f.write_all(contents.as_bytes())?;
    }

    let write_read = |out: &mut std::fs::File, r: &FastaRecord| -> anyhow::Result<()> {
        writeln!(out, ">{}\n{}", r.description, r.sequence)?;
        Ok(())
    };

    let mut next_written: Option<String> = None;
    for i in 0..records.len() {
        let curr = &records[i];
        let next = records.get(i + 1);
        if let Some(next) = next {
            let curr_parts: Vec<&str> = curr.id.splitn(2, '/').collect();
            let next_parts: Vec<&str> = next.id.splitn(2, '/').collect();
            if curr_parts.len() == 2
                && next_parts.len() == 2
                && curr_parts[0] == next_parts[0]
                && curr_parts[1] == "1"
                && next_parts[1] == "2"
            {
                write_read(&mut out_r1_f, curr)?;
                write_read(&mut out_r2_f, next)?;
                next_written = Some(next.id.clone());
            } else if Some(curr.id.clone()) != next_written {
                write_read(&mut out_rs_f, curr)?;
            }
        } else if Some(curr.id.clone()) != next_written {
            write_read(&mut out_rs_f, curr)?;
        }
    }
    Ok(())
}
