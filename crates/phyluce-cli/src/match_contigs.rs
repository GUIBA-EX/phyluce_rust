//! CLI wiring for `phyluce assembly match-contigs-to-probes`, mirroring
//! `phyluce_assembly_match_contigs_to_probes` end to end: LASTZ alignment of
//! each taxon's contigs against the probe set, duplicate contig/locus
//! filtering, and writing `probe.matches.sqlite` + optional CSV/dupe-report
//! output.

use std::collections::BTreeSet;
use std::io::Write;
use std::path::{Path, PathBuf};

use anyhow::Context;
use phyluce_assembly::{
    contig_count, contig_header_regex, contigs_matching_multiple_uces, db, extract_probe_name,
    get_probe_dupes, loci_matching_multiple_contigs, organism_names_from_fasta_paths,
    process_taxon_lastz_iter, FastMap, FastSet,
};
use phyluce_config::PhyluceConfig;
use phyluce_external::ExternalCommand;
use phyluce_io::lastz::{iter_lastz, read_lastz};
use phyluce_io::read_fasta;
use regex::Regex;

#[allow(clippy::too_many_arguments)]
pub fn run(
    contigs_dir: &Path,
    probes_path: &Path,
    output_dir: &Path,
    min_coverage: u32,
    min_identity: u32,
    dupefile: Option<PathBuf>,
    regex_str: &str,
    keep_duplicates: Option<PathBuf>,
    csv_path: Option<PathBuf>,
    skip_alignment: bool,
    force: bool,
) -> anyhow::Result<()> {
    // When --skip-alignment is set, the output dir is expected to already
    // hold precomputed `<stem>.lastz` files (e.g. replaying fixtures in an
    // environment without lastz installed), so it's left untouched.
    if skip_alignment {
        std::fs::create_dir_all(output_dir)
            .with_context(|| format!("creating output directory {}", output_dir.display()))?;
    } else if output_dir.exists() {
        if force {
            std::fs::remove_dir_all(output_dir)
                .with_context(|| format!("removing output directory {}", output_dir.display()))?;
            std::fs::create_dir_all(output_dir)
                .with_context(|| format!("creating output directory {}", output_dir.display()))?;
        } else {
            anyhow::bail!(
                "output directory {} already exists; pass --force to overwrite",
                output_dir.display()
            );
        }
    } else {
        std::fs::create_dir_all(output_dir)
            .with_context(|| format!("creating output directory {}", output_dir.display()))?;
    }

    let probe_regex = Regex::new(regex_str).context("invalid --regex")?;

    let probe_records = read_fasta(probes_path)
        .with_context(|| format!("reading probes fasta {}", probes_path.display()))?;
    let mut uces: BTreeSet<String> = BTreeSet::new();
    for record in &probe_records {
        uces.insert(extract_probe_name(&record.id, &probe_regex)?);
    }

    let dupes: FastSet<String> = match &dupefile {
        Some(path) => {
            let matches = read_lastz(path, false)
                .with_context(|| format!("reading dupefile {}", path.display()))?;
            get_probe_dupes(&matches, &probe_regex)?
        }
        None => FastSet::default(),
    };

    let mut fasta_files: Vec<PathBuf> = Vec::new();
    for ext in ["fasta", "fa", "fna"] {
        for entry in std::fs::read_dir(contigs_dir)
            .with_context(|| format!("reading contigs directory {}", contigs_dir.display()))?
        {
            let entry = entry?;
            let path = entry.path();
            if path.extension().and_then(|e| e.to_str()) == Some(ext) {
                fasta_files.push(path);
            }
        }
    }
    // Deterministic order for both file processing and (unlike the legacy
    // OS-dependent glob order) the SQLite column order -- see
    // docs/rust-rewrite-plan.md section 9's list of quirks worth fixing.
    fasta_files.sort();

    let organisms = organism_names_from_fasta_paths(&fasta_files).map_err(|errs| {
        let msgs: Vec<String> = errs.iter().map(|e| e.to_string()).collect();
        anyhow::anyhow!(msgs.join("\n"))
    })?;

    let cfg = PhyluceConfig::load()?;
    let header_fragments = cfg
        .get_contig_header_string()
        .context("no [headers] section in phyluce.conf")?;
    let header_regex = contig_header_regex(&header_fragments)?;

    let uces_vec: Vec<String> = uces.into_iter().collect();
    let db_path = output_dir.join("probe.matches.sqlite");
    if db_path.is_file() {
        std::fs::remove_file(&db_path)
            .with_context(|| format!("removing stale probe database {}", db_path.display()))?;
    }
    let conn = db::create_probe_database(&db_path, &organisms, &uces_vec)?;

    let mut dupe_writer = match &keep_duplicates {
        Some(p) => Some(
            std::fs::File::create(p)
                .with_context(|| format!("creating duplicates report file {}", p.display()))?,
        ),
        None => None,
    };
    let mut csv_writer = match &csv_path {
        Some(p) => {
            let mut w = csv::Writer::from_path(p)
                .with_context(|| format!("creating CSV output file {}", p.display()))?;
            w.write_record([
                "taxon",
                "uce-contigs",
                "total-contigs",
                "dupe-probe-matches",
                "loci-dropped",
                "contigs-dropped",
            ])?;
            Some(w)
        }
        None => None,
    };

    let lastz_bin = if skip_alignment {
        None
    } else {
        Some(cfg.get_user_path("binaries", "lastz")?)
    };

    for contig_path in &fasta_files {
        let file_stem = contig_path
            .file_stem()
            .and_then(|s| s.to_str())
            .unwrap_or("")
            .to_string();
        let critter = file_stem.split('.').next().unwrap_or("").replace('-', "_");
        let lastz_output = output_dir.join(format!(
            "{}.lastz",
            contig_path
                .file_stem()
                .and_then(|s| s.to_str())
                .unwrap_or("")
        ));
        let contigs = contig_count(contig_path)
            .with_context(|| format!("counting contigs in {}", contig_path.display()))?;

        if let Some(lastz_bin) = &lastz_bin {
            run_lastz_alignment(
                lastz_bin,
                contig_path,
                probes_path,
                min_coverage,
                min_identity,
                &lastz_output,
            )?;
        } else if !lastz_output.is_file() {
            anyhow::bail!(
                "--skip-alignment set but no pre-existing LASTZ output at {}",
                lastz_output.display()
            );
        }

        let result = process_taxon_lastz_iter(
            iter_lastz(&lastz_output, false)
                .with_context(|| format!("reading lastz output {}", lastz_output.display()))?,
            &probe_regex,
            &header_regex,
            &dupes,
            dupefile.is_some(),
        )?;

        if result.header_fallback_count > 0 {
            crate::cli_warn!(
                "{critter}: {} contig header(s) didn't match any [headers] pattern in \
                 phyluce.conf; used the header's first token as the contig name instead. \
                 Add a custom [headers] pattern if this looks wrong.",
                result.header_fallback_count
            );
        }

        let contigs_matching_mult_uces = contigs_matching_multiple_uces(&result.matches);
        let (uce_dupe_contigs, uce_dupe_uces) = loci_matching_multiple_contigs(&result.revmatches);
        let mut nodes_to_drop: FastSet<String> = contigs_matching_mult_uces.clone();
        nodes_to_drop.extend(uce_dupe_contigs.iter().cloned());

        if let Some(w) = dupe_writer.as_mut() {
            write_dupe_report(
                w,
                &critter,
                &uce_dupe_uces,
                &result.revmatches,
                &contigs_matching_mult_uces,
                &result.matches,
            )?;
        }

        let filtered: FastMap<String, FastSet<String>> = result
            .matches
            .iter()
            .filter(|(k, _)| !nodes_to_drop.contains(*k))
            .map(|(k, v)| (k.clone(), v.clone()))
            .collect();

        db::store_lastz_results(&conn, &filtered, &result.orientation, &critter)?;

        let unique_matches = filtered.len();
        let pct = if contigs > 0 {
            unique_matches as f64 / contigs as f64 * 100.0
        } else {
            0.0
        };
        crate::cli_info!(
            "{}: {} ({:.2}%) uniques of {} contigs, {} dupe probe matches, \
             {} UCE loci removed for matching multiple contigs, {} contigs \
             removed for matching multiple UCE loci",
            critter,
            unique_matches,
            pct,
            contigs,
            result.probe_dupes.len(),
            uce_dupe_uces.len(),
            contigs_matching_mult_uces.len(),
        );

        if let Some(w) = csv_writer.as_mut() {
            w.write_record([
                critter.as_str(),
                &unique_matches.to_string(),
                &contigs.to_string(),
                &result.probe_dupes.len().to_string(),
                &uce_dupe_uces.len().to_string(),
                &contigs_matching_mult_uces.len().to_string(),
            ])?;
        }
    }

    if let Some(mut w) = csv_writer {
        w.flush()?;
    }
    crate::cli_info!("The LASTZ alignments are in {}", output_dir.display());
    crate::cli_info!("The UCE match database is in {}", db_path.display());
    Ok(())
}

fn run_lastz_alignment(
    lastz_bin: &str,
    contig_path: &Path,
    probes_path: &Path,
    min_coverage: u32,
    min_identity: u32,
    output: &Path,
) -> anyhow::Result<()> {
    let target_arg = format!("{}[multiple,nameparse=full]", contig_path.display());
    let query_arg = format!("{}[nameparse=full]", probes_path.display());
    let report = ExternalCommand::new(lastz_bin)
        .args([
            target_arg,
            query_arg,
            "--strand=both".to_string(),
            "--seed=12of19".to_string(),
            "--transition".to_string(),
            "--nogfextend".to_string(),
            "--nochain".to_string(),
            "--gap=400,30".to_string(),
            "--xdrop=910".to_string(),
            "--ydrop=8370".to_string(),
            "--hspthresh=3000".to_string(),
            "--gappedthresh=3000".to_string(),
            "--noentropy".to_string(),
            format!("--coverage={min_coverage}"),
            format!("--identity={min_identity}"),
            format!("--output={}", output.display()),
            "--format=general-:score,name1,strand1,zstart1,end1,length1,name2,strand2,zstart2,end2,length2,diff,cigar,identity,continuity".to_string(),
        ])
        .run()
        .with_context(|| format!("running lastz for {}", contig_path.display()))?;
    if !report.stderr.trim().is_empty() {
        anyhow::bail!("lastz: {}", report.stderr);
    }
    Ok(())
}

fn write_dupe_report(
    w: &mut std::fs::File,
    critter: &str,
    uce_dupe_uces: &FastSet<String>,
    revmatches: &FastMap<String, FastSet<String>>,
    contigs_matching_mult_uces: &FastSet<String>,
    matches: &FastMap<String, FastSet<String>>,
) -> std::io::Result<()> {
    if !uce_dupe_uces.is_empty() {
        writeln!(w, "[{critter} - probes hitting multiple contigs]")?;
        let mut uces: Vec<&String> = uce_dupe_uces.iter().collect();
        uces.sort();
        for uce in uces {
            let mut contigs: Vec<&String> = revmatches.get(uce).into_iter().flatten().collect();
            contigs.sort();
            let joined = contigs
                .iter()
                .map(|s| s.as_str())
                .collect::<Vec<_>>()
                .join(", ");
            writeln!(w, "{uce}:{joined}")?;
        }
        writeln!(w)?;
    }
    if !contigs_matching_mult_uces.is_empty() {
        writeln!(w, "[{critter} - contigs hitting multiple probes]")?;
        let mut dupes: Vec<&String> = contigs_matching_mult_uces.iter().collect();
        dupes.sort();
        for dupe in dupes {
            let mut uces: Vec<&String> = matches.get(dupe).into_iter().flatten().collect();
            uces.sort();
            let joined = uces
                .iter()
                .map(|s| s.as_str())
                .collect::<Vec<_>>()
                .join(", ");
            writeln!(w, "{dupe}:{joined}")?;
        }
        writeln!(w)?;
    }
    Ok(())
}
