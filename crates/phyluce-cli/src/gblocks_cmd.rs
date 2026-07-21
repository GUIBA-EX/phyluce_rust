//! CLI wiring for `phyluce align get-gblocks-trimmed-alignments-from-untrimmed`,
//! mirroring `phyluce_align_get_gblocks_trimmed_alignments_from_untrimmed`.
//!
//! Gblocks only reads FASTA/NBRF-PIR, not NEXUS/PHYLIP/CLUSTAL/EMBOSS/
//! Stockholm -- true of the legacy Python script too (it hands Gblocks the
//! raw `--input-format` file unconverted, so `--input-format nexus` fails
//! there exactly the same way). Since the alignment is already parsed into
//! memory regardless of `--input-format`, this writes a temporary FASTA
//! copy for Gblocks to read instead of the original file, so every
//! `--input-format` this command accepts actually works with it.

use std::path::Path;

use anyhow::Context;
use phyluce_align::nexus::format_nexus;
use phyluce_config::PhyluceConfig;
use phyluce_io::read_fasta;

use crate::informative_sites_cmd::{find_alignment_files, load_alignment};

#[allow(clippy::too_many_arguments)]
pub fn run(
    alignments_dir: &Path,
    output_dir: &Path,
    input_format: &str,
    b1: f64,
    b2: f64,
    b3: u32,
    b4: u32,
    output_format: &str,
    cores: usize,
) -> anyhow::Result<()> {
    anyhow::ensure!(cores > 0, "--cores must be greater than zero");
    anyhow::ensure!(b4 >= 2, "--b4 must be >= 2");
    anyhow::ensure!(
        matches!(output_format, "fasta" | "nexus"),
        "output format '{output_format}' is not supported (only fasta/nexus)"
    );
    crate::output_path::prepare_output_dir(output_dir)?;

    let cfg = PhyluceConfig::load()?;
    let gblocks_bin = cfg.get_user_path("binaries", "gblocks")?;

    let files = find_alignment_files(alignments_dir, input_format)
        .with_context(|| format!("reading alignments directory {}", alignments_dir.display()))?;
    crate::parallel::ensure_unique_output_names(files.iter().map(|file| {
        let name = file
            .file_name()
            .and_then(|name| name.to_str())
            .unwrap_or("");
        let stem = name.split('.').next().unwrap_or(name);
        format!("{stem}.{output_format}")
    }))?;
    let count = files.len();
    let warnings = crate::parallel::try_map_ordered(files, cores, |file| {
        let name = file
            .file_name()
            .and_then(|n| n.to_str())
            .unwrap_or("")
            .split('.')
            .next()
            .unwrap_or("")
            .to_string();
        let alignment = load_alignment(&file, input_format)
            .with_context(|| format!("loading alignment {}", file.display()))?;
        let taxa = alignment.ntax();

        let b1_arg = (b1 * taxa as f64).round() as i64 + 1;
        let mut b2_arg = (b2 * taxa as f64).round() as i64;
        if b2_arg < b1_arg {
            b2_arg = b1_arg;
        }

        // Gblocks only reads FASTA/NBRF-PIR; write out the already-parsed
        // alignment as a temporary FASTA copy so `--input-format nexus`
        // (etc.) works instead of handing Gblocks a file format it can't
        // read (see the module doc comment).
        let gblocks_input = output_dir.join(format!(".{name}.gblocks-input.fasta"));
        // Gblocks writes beside the input using a fixed `-gb` suffix.
        let trimmed_path = std::path::PathBuf::from(format!("{}-gb", gblocks_input.display()));

        // Both temp files must be removed on every exit path -- success,
        // a Gblocks failure, or an error anywhere in between (a non-UTF8
        // sequence writing `gblocks_input`, `remove_stale_file` failing,
        // Gblocks output that fails to parse/validate, ...) -- not just
        // the specific paths that used to call `remove_file` individually.
        // Run the fallible work in a closure and clean up once afterward
        // instead of scattering removals along each early `?` return.
        let result = (|| -> anyhow::Result<Option<String>> {
            {
                let mut tmp = std::fs::File::create(&gblocks_input)
                    .with_context(|| format!("creating Gblocks FASTA input for locus {name}"))?;
                for row in &alignment.rows {
                    phyluce_io::write_fasta_record(
                        &mut tmp,
                        &row.id,
                        std::str::from_utf8(&row.seq)?,
                    )?;
                }
            }

            // Remove any prior result so a failed invocation cannot be
            // mistaken for a successful current run.
            crate::output_path::remove_stale_file(&trimmed_path)?;

            let output = std::process::Command::new(&gblocks_bin)
                .arg(&gblocks_input)
                .arg("-t=DNA")
                .arg(format!("-b1={b1_arg}"))
                .arg(format!("-b2={b2_arg}"))
                .arg(format!("-b3={b3}"))
                .arg(format!("-b4={b4}"))
                .arg("-b5=h")
                .arg("-p=n")
                .output()
                .with_context(|| format!("running Gblocks for locus {name}"))?;
            // Gblocks conventionally exits non-zero even on success; the
            // legacy script never checks the exit code, only whether the
            // `-gb` output file exists afterward.

            if !trimmed_path.is_file() {
                anyhow::ensure!(
                    output.status.success(),
                    "Gblocks failed for locus {name}: {}",
                    String::from_utf8_lossy(&output.stderr).trim()
                );
                return Ok(Some(format!("Missing information for locus {name}")));
            }
            let records = read_fasta(&trimmed_path)
                .with_context(|| format!("reading Gblocks output {}", trimmed_path.display()))?;
            let trimmed = phyluce_align::Alignment::from_pairs(
                records.into_iter().map(|r| (r.id, r.sequence)).collect(),
            );
            trimmed.validate()?;

            let ext = output_format;
            let out_path = output_dir.join(format!("{name}.{ext}"));
            if output_format == "fasta" {
                let mut out = std::fs::File::create(&out_path)
                    .with_context(|| format!("creating output file {}", out_path.display()))?;
                for row in &trimmed.rows {
                    phyluce_io::write_fasta_record(
                        &mut out,
                        &row.id,
                        std::str::from_utf8(&row.seq)?,
                    )?;
                }
            } else {
                std::fs::write(&out_path, format_nexus(&trimmed))
                    .with_context(|| format!("writing NEXUS output {}", out_path.display()))?;
            }
            Ok(None)
        })();
        let _ = std::fs::remove_file(&gblocks_input);
        let _ = std::fs::remove_file(&trimmed_path);
        result
    })?;
    for warning in warnings.into_iter().flatten() {
        crate::cli_warn!("{warning}");
    }
    for _ in 0..count {
        print!(".");
    }
    crate::cli_info!();
    Ok(())
}
