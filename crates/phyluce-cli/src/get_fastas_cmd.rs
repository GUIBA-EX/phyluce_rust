//! CLI wiring for `phyluce assembly get-fastas-from-match-counts`, mirroring
//! `phyluce_assembly_get_fastas_from_match_counts`.

use std::collections::HashSet;
use std::path::{Path, PathBuf};

use anyhow::Context;
use phyluce_assembly::get_fastas::{
    clean_sequence, find_contig_file, get_nodes_for_uces, node_with_strand_regex,
    reverse_complement,
};
use phyluce_assembly::match_counts::read_taxon_list_config;
use phyluce_assembly::{contig_header_regex, extract_contig_name};
use phyluce_config::PhyluceConfig;
use phyluce_io::{read_fasta, write_fasta_record};
use rusqlite::Connection;

#[allow(clippy::too_many_arguments)]
pub fn run(
    contigs_dir: &Path,
    locus_db: &Path,
    match_count_output: &Path,
    incomplete_matrix_out: Option<PathBuf>,
    output: &Path,
    extend_locus_db: Option<PathBuf>,
    extend_locus_contigs: Option<PathBuf>,
) -> anyhow::Result<()> {
    let config = read_taxon_list_config(match_count_output)?;
    let organisms = config
        .get("Organisms")
        .cloned()
        .context("no [Organisms] section in --match-count-output")?;
    let uces = config
        .get("Loci")
        .cloned()
        .context("no [Loci] section in --match-count-output")?;
    println!(
        "There are {} taxa in the match-count-config file named {}",
        organisms.len(),
        match_count_output.display()
    );
    if incomplete_matrix_out.is_some() {
        println!("There are {} UCE loci in an INCOMPLETE matrix", uces.len());
    } else {
        println!(
            "There are {} shared UCE loci in a COMPLETE matrix",
            uces.len()
        );
    }

    let conn = Connection::open(locus_db)?;
    if let Some(extend) = &extend_locus_db {
        conn.execute(
            &format!("ATTACH DATABASE '{}' AS extended", extend.display()),
            [],
        )?;
    }

    let cfg = PhyluceConfig::load()?;
    let header_fragments = cfg
        .get_contig_header_string()
        .context("no [headers] section in phyluce.conf")?;
    let header_regex = contig_header_regex(&header_fragments)?;
    let node_regex = node_with_strand_regex(&header_fragments)?;

    let mut incomplete_writer = match &incomplete_matrix_out {
        Some(p) => Some(std::fs::File::create(p)?),
        None => None,
    };

    let mut out = std::fs::File::create(output)?;
    let uces_all: HashSet<String> = uces.iter().cloned().collect();

    for organism in &organisms {
        println!("Getting UCE loci for {organism}");
        let name = organism.replace('_', "-");
        let is_extended = name.ends_with('*');

        let notstrict = incomplete_matrix_out.is_some();
        let (reads, node_dict, missing) = if !is_extended {
            let reads = find_contig_file(contigs_dir, &name)?;
            let (nd, missing) =
                get_nodes_for_uces(&conn, organism, &uces, false, notstrict, &node_regex)?;
            (reads, nd, missing)
        } else if let Some(extend_contigs) = &extend_locus_contigs {
            let stripped = name.trim_end_matches('*');
            let reads = find_contig_file(extend_contigs, stripped)?;
            let organism_stripped = organism.trim_end_matches('*');
            let (nd, missing) = get_nodes_for_uces(
                &conn,
                organism_stripped,
                &uces,
                true,
                notstrict,
                &node_regex,
            )?;
            (reads, nd, missing)
        } else {
            anyhow::bail!(
                "organism '{organism}' is extended (trailing '*') but --extend-locus-contigs was not given"
            );
        };

        println!("There are {} UCE loci for {organism}", node_dict.len());
        println!("Parsing and renaming contigs for {organism}");

        let node_dict_set: HashSet<&String> = node_dict.keys().collect();
        let mut nodes_written: HashSet<String> = HashSet::new();
        let mut written: Vec<String> = Vec::new();
        let mut n_replace_count = 0usize;

        let records = read_fasta(&reads)?;
        for record in &records {
            if nodes_written.len() == node_dict_set.len() {
                break;
            }
            let contig_name = extract_contig_name(&record.id, &header_regex)?.to_lowercase();
            if let Some(node) = node_dict.get(&contig_name) {
                let organism_stripped = organism.trim_end_matches('*');
                let new_id = format!("{}_{} |{}", node.uce, organism_stripped, node.uce);
                let (mut cleaned, replaced) = clean_sequence(&record.sequence);
                if replaced {
                    n_replace_count += 1;
                }
                if node.strand == '-' {
                    cleaned = reverse_complement(&cleaned);
                }
                write_fasta_record(&mut out, &new_id, &cleaned)?;
                written.push(node.uce.clone());
                nodes_written.insert(contig_name);
            }
        }

        if n_replace_count > 0 {
            println!(
                "Replaced <20 ambiguous bases (N) in {n_replace_count} contigs for {organism}"
            );
        }

        if let (Some(w), false) = (incomplete_writer.as_mut(), missing.is_empty()) {
            println!(
                "Writing missing locus information to {}",
                incomplete_matrix_out.as_ref().unwrap().display()
            );
            use std::io::Write as _;
            writeln!(w, "[{organism}]")?;
            for m in &missing {
                writeln!(w, "{m}")?;
                written.push(m.clone());
            }
        }

        let written_set: HashSet<String> = written.into_iter().collect();
        anyhow::ensure!(
            written_set == uces_all,
            "UCE names do not match for organism {organism}"
        );
    }
    Ok(())
}
