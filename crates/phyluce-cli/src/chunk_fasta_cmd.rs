//! CLI wiring for `phyluce ncbi chunk-fasta-for-ncbi`, mirroring
//! `phyluce_ncbi_chunk_fasta_for_ncbi`.

use std::path::Path;

use phyluce_io::{read_fasta, write_fasta_record};

pub fn run(
    input: &Path,
    chunk_size: usize,
    output_prefix: &str,
    output_suffix: &str,
) -> anyhow::Result<()> {
    anyhow::ensure!(chunk_size > 0, "--chunk-size must be greater than zero");
    let records = read_fasta(input)?;
    for (i, batch) in records.chunks(chunk_size).enumerate() {
        let filename = format!("{output_prefix}_{}.{output_suffix}", i + 1);
        let mut out = std::fs::File::create(&filename)?;
        for record in batch {
            write_fasta_record(&mut out, &record.description, &record.sequence)?;
        }
        crate::cli_info!("Wrote {} records to {filename}", batch.len());
    }
    Ok(())
}
