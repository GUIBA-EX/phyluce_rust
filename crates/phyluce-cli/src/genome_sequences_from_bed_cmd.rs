//! CLI wiring for `phyluce probe get-genome-sequences-from-bed`, mirroring
//! `phyluce_probe_get_genome_sequences_from_bed`.

use std::io::Write as _;
use std::path::Path;

use anyhow::Context;
use phyluce_io::twobit::TwoBitFile;

pub fn run(
    bed: &Path,
    twobit: &Path,
    output: &Path,
    filter_mask: Option<f64>,
    max_n: usize,
    buffer_to: Option<i64>,
) -> anyhow::Result<()> {
    let tb = TwoBitFile::open(twobit)
        .with_context(|| format!("opening 2bit file {}", twobit.display()))?;
    let bed_text = std::fs::read_to_string(bed)
        .with_context(|| format!("reading BED file {}", bed.display()))?;
    let mut out = std::fs::File::create(output)
        .with_context(|| format!("creating output file {}", output.display()))?;

    let mut filtered = 0usize;
    let mut kept = 0usize;
    let mut total_lines = 0usize;
    for (cnt, line) in bed_text.lines().enumerate() {
        total_lines += 1;
        let fields: Vec<&str> = line.trim().split('\t').collect();
        anyhow::ensure!(fields.len() == 3, "malformed BED line: {line:?}");
        let chromo = fields[0];
        let mut start: i64 = fields[1].parse()?;
        let mut end: i64 = fields[2].parse()?;
        if let Some(buffer_to) = buffer_to {
            let length = (end - start).abs();
            let mut delta = buffer_to - length;
            if delta > 0 {
                if delta % 2 != 0 {
                    delta += 1;
                }
                start -= delta / 2;
                end += delta / 2;
            }
        }
        let sequence = tb.read_slice(chromo, start, end)?;
        let n_count = sequence
            .iter()
            .filter(|&&b| b.eq_ignore_ascii_case(&b'N'))
            .count();
        let masked = sequence.iter().filter(|b| b.is_ascii_lowercase()).count() as f64
            / sequence.len().max(1) as f64;
        let is_masked = filter_mask.map(|m| masked > m).unwrap_or(false);
        let long_enough = buffer_to
            .map(|b| sequence.len() as i64 >= b)
            .unwrap_or(true);
        if long_enough && n_count <= max_n && !is_masked {
            writeln!(
                out,
                ">slice_{cnt} |{chromo}:{start}-{end}\n{}",
                String::from_utf8_lossy(&sequence)
            )?;
            kept += 1;
        } else {
            filtered += 1;
        }
    }

    // Not `cnt + 1` (Python's approach, using the *last* loop index):
    // that only works because Python's `enumerate()` also happens to
    // crash with `UnboundLocalError` on an empty BED file (`cnt` is never
    // bound), so it never has to handle `cnt == 0` meaning zero lines
    // instead of one. Track the real line count instead, so an empty BED
    // file cleanly reports 0 rather than misreporting 1.
    crate::cli_warn!(
        "Screened {total_lines} sequences.  Filtered {filtered} < 160 bp or with > {}% masked bases or > {max_n} N-bases. Kept {kept}.",
        filter_mask.unwrap_or(0.25) * 100.0,
    );
    Ok(())
}
