//! CLI wiring for `phyluce assembly match-contigs-to-barcodes`, mirroring
//! `phyluce_assembly_match_contigs_to_barcodes`.
//!
//! Untested: needs the `lastz` binary (not installed here). Also, the
//! BOLD systems web-lookup step (`requests.get("http://v4.boldsystems.org/...")`)
//! is NOT reproduced. Callers must pass `--no-bold`; otherwise this command
//! errors instead of silently skipping the lookup.

use std::collections::HashMap;
use std::path::Path;

use phyluce_config::PhyluceConfig;
use phyluce_external::ExternalCommand;
use phyluce_io::{read_fasta, write_fasta_record, FastaRecord};

fn reverse_complement(seq: &str) -> String {
    seq.chars()
        .rev()
        .map(|c| match c.to_ascii_uppercase() {
            'A' => 'T',
            'T' => 'A',
            'C' => 'G',
            'G' => 'C',
            other => other,
        })
        .collect()
}

pub fn run(
    contigs_dir: &Path,
    barcodes: &Path,
    output_dir: &Path,
    no_bold: bool,
    _database: &str,
) -> anyhow::Result<()> {
    anyhow::ensure!(
        no_bold,
        "BOLD lookups are not implemented in the Rust port; pass --no-bold to run LASTZ slicing only"
    );
    std::fs::create_dir_all(output_dir)?;
    let cfg = PhyluceConfig::load()?;
    let lastz_bin = cfg.get_user_path("binaries", "lastz")?;

    let mut fasta_files: Vec<std::path::PathBuf> = std::fs::read_dir(contigs_dir)?
        .filter_map(|e| e.ok())
        .map(|e| e.path())
        .filter(|p| {
            p.extension()
                .and_then(|e| e.to_str())
                .is_some_and(|e| e.starts_with("fa"))
        })
        .collect();
    fasta_files.sort();

    for contig_path in &fasta_files {
        let records = read_fasta(contig_path)?;
        let record_dict: HashMap<String, FastaRecord> =
            records.into_iter().map(|r| (r.id.clone(), r)).collect();

        let stem = contig_path
            .file_stem()
            .and_then(|s| s.to_str())
            .unwrap_or("");
        let lastz_output = output_dir.join(format!("{stem}.lastz"));

        ExternalCommand::new(&lastz_bin)
            .args([
                format!("{}[multiple,nameparse=full]", barcodes.display()),
                format!("{}[nameparse=full]", contig_path.display()),
                format!("--output={}", lastz_output.display()),
                "--format=general-:score,name1,strand1,zstart1,end1,length1,name2,strand2,zstart2,end2,length2,diff,cigar,identity,continuity".to_string(),
            ])
            .run()?;

        let matches = phyluce_io::lastz::read_lastz(&lastz_output, false)?;
        let mut slices: Vec<(String, String)> = Vec::new();
        for m in &matches {
            let matching_contig = m.name2.split(' ').next().unwrap_or(&m.name2).to_string();
            let Some(record) = record_dict.get(&matching_contig) else {
                eprintln!("Did not find a match for locus {matching_contig}");
                continue;
            };
            let start = m.zstart2 as usize;
            let end = m.end2 as usize;
            let seq = if m.strand2 == "+" {
                record.sequence.get(start..end).unwrap_or("").to_string()
            } else {
                let rc = reverse_complement(&record.sequence);
                rc.get(start..end).unwrap_or("").to_string()
            };
            slices.push((matching_contig, seq));
        }

        let slices_output = output_dir.join(format!("{stem}.slices.fasta"));
        let mut outf = std::fs::File::create(&slices_output)?;
        for (id, seq) in &slices {
            write_fasta_record(&mut outf, id, seq)?;
        }
    }
    Ok(())
}
