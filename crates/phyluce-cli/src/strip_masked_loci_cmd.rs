//! CLI wiring for `phyluce probe strip-masked-loci-from-set`, mirroring
//! `phyluce_probe_strip_masked_loci_from_set`.

use std::io::Write as _;
use std::path::Path;

use phyluce_io::twobit::TwoBitFile;

pub fn run(
    bed: &Path,
    twobit: &Path,
    output: &Path,
    filter_mask: Option<f64>,
    max_n: usize,
    min_length: i64,
) -> anyhow::Result<()> {
    let tb = TwoBitFile::open(twobit)?;
    let bed_text = std::fs::read_to_string(bed)?;
    let mut out = std::fs::File::create(output)?;

    let mut filtered = 0usize;
    let mut kept = 0usize;
    let mut cnt = 0usize;
    for (i, line) in bed_text.lines().enumerate() {
        cnt = i;
        let fields: Vec<&str> = line.trim().split('\t').collect();
        anyhow::ensure!(fields.len() == 3, "malformed BED line: {line:?}");
        let chromo = fields[0];
        let start: i64 = fields[1].parse()?;
        let end: i64 = fields[2].parse()?;
        let sequence = tb.read_slice(chromo, start, end)?;
        let n_count = sequence
            .iter()
            .filter(|&&b| b.eq_ignore_ascii_case(&b'N'))
            .count();
        let masked = sequence.iter().filter(|b| b.is_ascii_lowercase()).count() as f64
            / sequence.len().max(1) as f64;
        let is_masked = filter_mask.map(|m| masked > m).unwrap_or(false);
        let long_enough = min_length <= 0 || (end - start) >= min_length;
        if n_count <= max_n && !is_masked && long_enough {
            writeln!(out, "{chromo}\t{start}\t{end}")?;
            kept += 1;
        } else {
            filtered += 1;
        }
    }

    crate::cli_warn!(
        "Screened {} sequences from {}.  Filtered {filtered} with > {}% masked bases or > {max_n} N-bases or < {min_length} length. Kept {kept}.",
        cnt + 1,
        bed.file_name().and_then(|s| s.to_str()).unwrap_or(""),
        filter_mask.unwrap_or(0.0) * 100.0,
    );
    Ok(())
}
