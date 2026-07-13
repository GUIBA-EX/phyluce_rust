//! CLI wiring for `phyluce assembly explode-get-fastas-file`, mirroring
//! `phyluce_assembly_explode_get_fastas_file`.

use std::collections::HashMap;
use std::path::{Path, PathBuf};

use phyluce_assembly::explode::{locus_key, taxon_key};
use phyluce_io::{read_fasta, write_fasta_record};

pub fn run(input: &Path, output: &Path, by_taxon: bool, split_char: &str) -> anyhow::Result<()> {
    std::fs::create_dir_all(output)?;

    crate::cli_info!("Reading fasta...");
    let records = read_fasta(input)?;
    let mut groups: HashMap<String, Vec<&phyluce_io::FastaRecord>> = HashMap::new();
    for record in &records {
        let key = if by_taxon {
            taxon_key(&record.id, split_char)
        } else {
            locus_key(&record.id, split_char)
        };
        groups.entry(key).or_default().push(record);
    }

    crate::cli_info!("Writing fasta...");
    for (key, seqs) in &groups {
        let out_path: PathBuf =
            crate::output_path::output_file(output, &format!("{key}.unaligned.fasta"))?;
        let mut out = std::fs::File::create(out_path)?;
        for seq in seqs {
            write_fasta_record(&mut out, &seq.description, &seq.sequence)?;
        }
    }
    Ok(())
}
