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
use phyluce_assembly::{contig_header_regex, extract_contig_name_lenient};
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
    let config = read_taxon_list_config(match_count_output).with_context(|| {
        format!(
            "reading match-count-output config {}",
            match_count_output.display()
        )
    })?;
    let organisms = config
        .get("Organisms")
        .cloned()
        .context("no [Organisms] section in --match-count-output")?;
    let uces = config
        .get("Loci")
        .cloned()
        .context("no [Loci] section in --match-count-output")?;
    crate::cli_info!(
        "There are {} taxa in the match-count-config file named {}",
        organisms.len(),
        match_count_output.display()
    );
    if incomplete_matrix_out.is_some() {
        crate::cli_info!("There are {} UCE loci in an INCOMPLETE matrix", uces.len());
    } else {
        crate::cli_info!(
            "There are {} shared UCE loci in a COMPLETE matrix",
            uces.len()
        );
    }

    let conn = Connection::open(locus_db)
        .with_context(|| format!("opening locus database {}", locus_db.display()))?;
    if let Some(extend) = &extend_locus_db {
        conn.execute(
            "ATTACH DATABASE ?1 AS extended",
            [extend.to_string_lossy().as_ref()],
        )
        .with_context(|| format!("attaching extended locus database {}", extend.display()))?;
    }

    let cfg = PhyluceConfig::load()?;
    let header_fragments = cfg
        .get_contig_header_string()
        .context("no [headers] section in phyluce.conf")?;
    let header_regex = contig_header_regex(&header_fragments)?;
    let node_regex = node_with_strand_regex(&header_fragments)?;

    let mut incomplete_writer =
        match &incomplete_matrix_out {
            Some(p) => Some(std::fs::File::create(p).with_context(|| {
                format!("creating incomplete matrix output file {}", p.display())
            })?),
            None => None,
        };

    let mut out = std::fs::File::create(output)
        .with_context(|| format!("creating output file {}", output.display()))?;
    let uces_all: HashSet<String> = uces.iter().cloned().collect();

    for organism in &organisms {
        crate::cli_info!("Getting UCE loci for {organism}");
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

        crate::cli_info!("There are {} UCE loci for {organism}", node_dict.len());
        crate::cli_info!("Parsing and renaming contigs for {organism}");

        let node_dict_set: HashSet<&String> = node_dict.keys().collect();
        let mut nodes_written: HashSet<String> = HashSet::new();
        let mut written: Vec<String> = Vec::new();
        let mut n_replace_count = 0usize;
        let mut header_fallback_count = 0usize;

        let records = read_fasta(&reads)
            .with_context(|| format!("reading contigs fasta {}", reads.display()))?;
        for record in &records {
            if nodes_written.len() == node_dict_set.len() {
                break;
            }
            let (contig_name, used_fallback) =
                extract_contig_name_lenient(&record.id, &header_regex)?;
            let contig_name = contig_name.to_lowercase();
            if used_fallback {
                header_fallback_count += 1;
            }
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
            crate::cli_info!(
                "Replaced <20 ambiguous bases (N) in {n_replace_count} contigs for {organism}"
            );
        }
        if header_fallback_count > 0 {
            crate::cli_warn!(
                "{organism}: {header_fallback_count} contig header(s) didn't match any \
                 [headers] pattern in phyluce.conf; used the header's first token as the \
                 contig name instead. Add a custom [headers] pattern if this looks wrong."
            );
        }

        if let (Some(w), false) = (incomplete_writer.as_mut(), missing.is_empty()) {
            crate::cli_info!(
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
