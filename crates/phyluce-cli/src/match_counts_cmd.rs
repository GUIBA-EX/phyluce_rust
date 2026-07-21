//! CLI wiring for `phyluce assembly get-match-counts`.

use std::path::{Path, PathBuf};
use std::time::{SystemTime, UNIX_EPOCH};

use anyhow::Context;
use phyluce_assembly::match_counts::{
    complete_matrix, format_output, incomplete_matrix, matches_by_organism, read_taxon_list_config,
    sample_optimized_groups, taxa_from_config, uce_names, OptimizedGroup,
};
use rusqlite::Connection;

#[derive(Clone, Copy, Debug)]
pub struct OptimizationOptions {
    pub optimize: bool,
    pub random: bool,
    pub samples: usize,
    pub sample_size: usize,
    pub silent: bool,
    pub keep_counts: bool,
    pub seed: Option<u64>,
    pub cores: usize,
}

fn generated_seed() -> u64 {
    let nanos = SystemTime::now()
        .duration_since(UNIX_EPOCH)
        .unwrap_or_default()
        .as_nanos();
    nanos as u64 ^ (nanos >> 64) as u64 ^ (u64::from(std::process::id()) << 32)
}

fn format_exhaustive_report(results: &[OptimizedGroup]) -> String {
    let mut report = results
        .iter()
        .map(|result| format!("{}\t{}", result.organisms.join(","), result.locus_count()))
        .collect::<Vec<_>>()
        .join("\n");
    report.push('\n');
    report
}

fn run_exhaustive_optimization(
    organismal_matches: &std::collections::HashMap<String, std::collections::HashSet<String>>,
    organisms: &[String],
    uces: &std::collections::HashSet<String>,
    output: &Path,
    cores: usize,
) -> anyhow::Result<()> {
    anyhow::ensure!(
        !organisms.is_empty(),
        "cannot optimize an empty taxon group"
    );
    anyhow::ensure!(cores > 0, "--cores must be greater than zero");
    crate::cli_info!(
        "Enumerating the best complete-matrix group for sizes 1 through {} using {} worker(s)",
        organisms.len(),
        cores.min(organisms.len())
    );
    let sizes: Vec<usize> = (1..=organisms.len()).collect();
    let results = crate::parallel::try_map_ordered(sizes, cores, |size| {
        Ok(phyluce_assembly::match_counts::optimize_group_for_size(
            organismal_matches,
            organisms,
            uces,
            size,
        )?)
    })?;
    let report = format_exhaustive_report(&results);
    print!("{report}");
    tracing::info!(message = %report.trim_end());
    crate::output_path::write_atomic(output, report)?;
    Ok(())
}

fn run_random_optimization(
    organismal_matches: &std::collections::HashMap<String, std::collections::HashSet<String>>,
    organisms: &[String],
    uces: &std::collections::HashSet<String>,
    output: &Path,
    options: OptimizationOptions,
) -> anyhow::Result<()> {
    let seed = options.seed.unwrap_or_else(generated_seed);
    crate::cli_info!(
        "Sampling {} candidate groups of {} taxa with random seed {seed}",
        options.samples,
        options.sample_size
    );
    let result = sample_optimized_groups(
        organismal_matches,
        organisms,
        uces,
        options.samples,
        options.sample_size,
        seed,
    )?;

    if options.keep_counts {
        let counts = result
            .counts
            .iter()
            .map(|(size, count)| format!("{size},{count}"))
            .collect::<Vec<_>>()
            .join("\n");
        crate::output_path::write_atomic(output, counts)?;
        return Ok(());
    }

    let mut best_group = result.best.organisms.clone();
    best_group.sort();
    let missing = result
        .missing_counts
        .iter()
        .map(|(taxon, count)| format!("{taxon}: {count}"))
        .collect::<Vec<_>>()
        .join(", ");
    crate::cli_info!("max UCE = {}", result.best.locus_count());
    crate::cli_info!("group size = {}", result.best.organisms.len());
    crate::cli_info!("best group\n\t{:?}", best_group);
    crate::cli_info!("Times not in best group per iteration\n\t{{{missing}}}");

    if !options.silent {
        crate::cli_info!(
            "Writing the optimized taxa and loci in the data matrix to {}",
            output.display()
        );
        crate::output_path::write_atomic(
            output,
            format_output(&result.best.organisms, &result.best.uces),
        )?;
    }
    Ok(())
}

#[allow(clippy::too_many_arguments)]
pub fn run(
    locus_db: &Path,
    taxon_list_config: &Path,
    taxon_group: &str,
    output: &Path,
    incomplete_matrix_flag: bool,
    extend_locus_db: Option<PathBuf>,
    options: OptimizationOptions,
) -> anyhow::Result<()> {
    anyhow::ensure!(
        !options.random || options.optimize,
        "--random requires --optimize"
    );
    anyhow::ensure!(
        !options.keep_counts || (options.optimize && options.random),
        "--keep-counts requires --optimize --random"
    );
    anyhow::ensure!(
        !options.optimize || !incomplete_matrix_flag,
        "--optimize maximizes complete-matrix loci and cannot be combined with --incomplete-matrix"
    );
    let mut protected_inputs = vec![locus_db, taxon_list_config];
    if let Some(extend) = &extend_locus_db {
        protected_inputs.push(extend);
    }
    crate::output_path::ensure_output_not_input(output, &protected_inputs)?;

    let conn = Connection::open(locus_db)
        .with_context(|| format!("opening locus database {}", locus_db.display()))?;
    if let Some(extend) = &extend_locus_db {
        conn.execute(
            "ATTACH DATABASE ?1 AS extended",
            [extend.to_string_lossy().as_ref()],
        )?;
    }

    let config = read_taxon_list_config(taxon_list_config)
        .with_context(|| format!("reading taxon list config {}", taxon_list_config.display()))?;
    let organisms = taxa_from_config(&config, taxon_group)
        .with_context(|| format!("taxon-group '{taxon_group}'"))?;
    crate::cli_info!(
        "There are {} taxa in the taxon-group '[{}]' in the config file {}",
        organisms.len(),
        taxon_group,
        taxon_list_config.display()
    );

    let uces = uce_names(&conn)?;
    crate::cli_info!("There are {} total UCE loci in the database", uces.len());

    let organismal_matches = matches_by_organism(&conn, &organisms)?;
    if options.optimize {
        if options.random {
            return run_random_optimization(
                &organismal_matches,
                &organisms,
                &uces,
                output,
                options,
            );
        }
        return run_exhaustive_optimization(
            &organismal_matches,
            &organisms,
            &uces,
            output,
            options.cores,
        );
    }

    let shared_uces = if !incomplete_matrix_flag {
        crate::cli_info!("Getting UCE matches by organism to generate a COMPLETE matrix");
        let (shared, losses) = complete_matrix(&organismal_matches, &organisms, &uces);
        crate::cli_info!(
            "There are {} shared UCE loci in a COMPLETE matrix",
            shared.len()
        );
        let mut sorted_losses: Vec<(&String, &usize)> = losses.iter().collect();
        sorted_losses.sort_by(|a, b| b.1.cmp(a.1));
        for (organism, loss) in sorted_losses {
            crate::cli_info!("\tFailed to detect {loss} UCE loci in {organism}");
        }
        shared
    } else {
        crate::cli_info!("Getting UCE matches by organism to generate a INCOMPLETE matrix");
        let shared = incomplete_matrix(&organismal_matches, &uces);
        crate::cli_info!(
            "There are {} UCE loci in an INCOMPLETE matrix",
            shared.len()
        );
        shared
    };

    if !organisms.is_empty() && !options.silent {
        crate::cli_info!(
            "Writing the taxa and loci in the data matrix to {}",
            output.display()
        );
        crate::output_path::write_atomic(output, format_output(&organisms, &shared_uces))?;
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::sync::atomic::{AtomicUsize, Ordering};

    static TEST_COUNTER: AtomicUsize = AtomicUsize::new(0);

    fn fixture() -> (PathBuf, PathBuf, PathBuf) {
        let root = std::env::temp_dir().join(format!(
            "phyluce-match-count-optimize-{}-{}",
            std::process::id(),
            TEST_COUNTER.fetch_add(1, Ordering::Relaxed)
        ));
        std::fs::create_dir_all(&root).unwrap();
        let database = root.join("matches.sqlite");
        let config = root.join("taxa.conf");
        let output = root.join("output.conf");
        let conn = Connection::open(&database).unwrap();
        conn.execute_batch(
            "
            CREATE TABLE matches (uce TEXT PRIMARY KEY, a INTEGER, b INTEGER, c INTEGER);
            CREATE TABLE match_map (uce TEXT PRIMARY KEY, a TEXT, b TEXT, c TEXT);
            INSERT INTO matches VALUES
                ('uce-1', 1, 1, 1),
                ('uce-2', 1, 1, 0),
                ('uce-3', 1, 1, 0),
                ('uce-4', 1, 0, 1);
            INSERT INTO match_map VALUES
                ('uce-1', 'a1(+)', 'b1(+)', 'c1(+)'),
                ('uce-2', 'a2(+)', 'b2(+)', ''),
                ('uce-3', 'a3(+)', 'b3(+)', ''),
                ('uce-4', 'a4(+)', '', 'c4(+)');
            ",
        )
        .unwrap();
        drop(conn);
        std::fs::write(&config, "[all]\na\nb\nc\n").unwrap();
        (database, config, output)
    }

    fn options() -> OptimizationOptions {
        OptimizationOptions {
            optimize: true,
            random: false,
            samples: 10,
            sample_size: 2,
            silent: false,
            keep_counts: false,
            seed: Some(7),
            cores: 2,
        }
    }

    #[test]
    fn exhaustive_optimization_writes_each_best_group() {
        let (database, config, output) = fixture();
        run(&database, &config, "all", &output, false, None, options()).unwrap();
        assert_eq!(
            std::fs::read_to_string(&output).unwrap(),
            "a\t4\na,b\t3\na,b,c\t1\n"
        );
        let serial_output = output.with_file_name("serial.tsv");
        let mut serial_options = options();
        serial_options.cores = 1;
        run(
            &database,
            &config,
            "all",
            &serial_output,
            false,
            None,
            serial_options,
        )
        .unwrap();
        assert_eq!(
            std::fs::read_to_string(&serial_output).unwrap(),
            std::fs::read_to_string(&output).unwrap()
        );
        std::fs::remove_dir_all(database.parent().unwrap()).unwrap();
    }

    #[test]
    fn random_optimization_writes_best_matrix_with_a_seed() {
        let (database, config, output) = fixture();
        let mut options = options();
        options.random = true;
        options.samples = 3;
        run(&database, &config, "all", &output, false, None, options).unwrap();
        assert_eq!(
            std::fs::read_to_string(&output).unwrap(),
            "[Organisms]\na\nb\n[Loci]\nuce-1\nuce-2\nuce-3\n"
        );
        std::fs::remove_dir_all(database.parent().unwrap()).unwrap();
    }

    #[test]
    fn random_keep_counts_writes_one_row_per_iteration() {
        let (database, config, output) = fixture();
        let mut options = options();
        options.random = true;
        options.keep_counts = true;
        options.samples = 3;
        run(&database, &config, "all", &output, false, None, options).unwrap();
        assert_eq!(std::fs::read_to_string(&output).unwrap(), "2,3\n2,3\n2,3");
        std::fs::remove_dir_all(database.parent().unwrap()).unwrap();
    }

    #[test]
    fn refuses_to_overwrite_the_locus_database() {
        let (database, config, _) = fixture();
        let error = run(&database, &config, "all", &database, false, None, options()).unwrap_err();
        assert!(error.to_string().contains("must not overwrite input"));
        assert_eq!(
            Connection::open(&database)
                .unwrap()
                .query_row("SELECT COUNT(*) FROM matches", [], |row| row
                    .get::<_, i64>(0))
                .unwrap(),
            4
        );
        std::fs::remove_dir_all(database.parent().unwrap()).unwrap();
    }
}
