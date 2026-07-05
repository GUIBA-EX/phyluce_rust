//! CLI wiring for `phyluce assembly get-match-counts`, mirroring
//! `phyluce_assembly_get_match_counts`'s non-`--optimize` path.

use std::path::{Path, PathBuf};

use anyhow::Context;
use phyluce_assembly::match_counts::{
    complete_matrix, format_output, incomplete_matrix, matches_by_organism, read_taxon_list_config,
    taxa_from_config, uce_names,
};
use rusqlite::Connection;

#[allow(clippy::too_many_arguments)]
pub fn run(
    locus_db: &Path,
    taxon_list_config: &Path,
    taxon_group: &str,
    output: &Path,
    incomplete_matrix_flag: bool,
    extend_locus_db: Option<PathBuf>,
) -> anyhow::Result<()> {
    let conn = Connection::open(locus_db)?;
    if let Some(extend) = &extend_locus_db {
        conn.execute(
            &format!("ATTACH DATABASE '{}' AS extended", extend.display()),
            [],
        )?;
    }

    let config = read_taxon_list_config(taxon_list_config)?;
    let organisms = taxa_from_config(&config, taxon_group)
        .with_context(|| format!("taxon-group '{taxon_group}'"))?;
    println!(
        "There are {} taxa in the taxon-group '[{}]' in the config file {}",
        organisms.len(),
        taxon_group,
        taxon_list_config.display()
    );

    let uces = uce_names(&conn)?;
    println!("There are {} total UCE loci in the database", uces.len());

    let organismal_matches = matches_by_organism(&conn, &organisms)?;
    let shared_uces = if !incomplete_matrix_flag {
        println!("Getting UCE matches by organism to generate a COMPLETE matrix");
        let (shared, losses) = complete_matrix(&organismal_matches, &organisms, &uces);
        println!(
            "There are {} shared UCE loci in a COMPLETE matrix",
            shared.len()
        );
        let mut sorted_losses: Vec<(&String, &usize)> = losses.iter().collect();
        sorted_losses.sort_by(|a, b| b.1.cmp(a.1));
        for (organism, loss) in sorted_losses {
            println!("\tFailed to detect {loss} UCE loci in {organism}");
        }
        shared
    } else {
        println!("Getting UCE matches by organism to generate a INCOMPLETE matrix");
        let shared = incomplete_matrix(&organismal_matches, &uces);
        println!(
            "There are {} UCE loci in an INCOMPLETE matrix",
            shared.len()
        );
        shared
    };

    if !organisms.is_empty() {
        println!(
            "Writing the taxa and loci in the data matrix to {}",
            output.display()
        );
        std::fs::write(output, format_output(&organisms, &shared_uces))?;
    }
    Ok(())
}
