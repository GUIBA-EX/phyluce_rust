//! CLI wiring for `phyluce assembly extract-contigs-to-barcodes`, mirroring
//! `phyluce_assembly_extract_contigs_to_barcodes`.

use std::collections::HashMap;
use std::path::Path;

use phyluce_io::{read_fasta, write_fasta_record, FastaRecord};

pub fn run(contigs_dir: &Path, config: &Path, output: &Path) -> anyhow::Result<()> {
    let config_text = std::fs::read_to_string(config)?;
    let ini = phyluce_config::Ini::parse(&config_text);
    let entries = ini
        .entries("assemblies")
        .ok_or_else(|| anyhow::anyhow!("no [assemblies] section in --config"))?;

    let mut outf = std::fs::File::create(output)?;
    let mut cache: HashMap<String, Option<HashMap<String, FastaRecord>>> = HashMap::new();

    for (raw_key, contig) in entries {
        let assembly = raw_key.split('|').next().unwrap_or(raw_key).to_string();
        let critter = assembly
            .split('.')
            .next()
            .unwrap_or(&assembly)
            .replace('-', "_");

        if !cache.contains_key(&assembly) {
            let path = contigs_dir.join(assembly.replace('_', "-"));
            let index = match read_fasta(&path) {
                Ok(records) => Some(records.into_iter().map(|r| (r.id.clone(), r)).collect()),
                Err(_) => None,
            };
            cache.insert(assembly.clone(), index);
        }

        if let Some(Some(index)) = cache.get(&assembly) {
            let record = index.get(contig).ok_or_else(|| {
                anyhow::anyhow!("contig '{contig}' not found in assembly '{assembly}'")
            })?;
            write_fasta_record(
                &mut outf,
                &format!("{critter}|{}", record.id),
                &record.sequence,
            )?;
        }
    }
    Ok(())
}
