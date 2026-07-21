//! CLI wiring for `phyluce align concatenate-alignments`, mirroring
//! `phyluce_align_concatenate_alignments` with a disk-backed matrix.

use std::collections::{HashMap, HashSet};
use std::fs::{File, OpenOptions};
use std::io::{BufWriter, Read, Seek, SeekFrom, Write};
use std::path::{Path, PathBuf};
use std::sync::atomic::{AtomicU64, Ordering};

use anyhow::Context;
use phyluce_align::concat::{format_sets_block, Charset};
use phyluce_align::nexus::safename;

use crate::informative_sites_cmd::{find_alignment_files, load_alignment};

static TEMP_COUNTER: AtomicU64 = AtomicU64::new(0);

struct ScratchMatrix {
    path: PathBuf,
    file: File,
}

struct ConcatenationMetadata {
    taxa: Vec<String>,
    taxon_index: HashMap<String, usize>,
    charsets: Vec<Charset>,
    locus_taxa: Vec<HashSet<String>>,
    total_length: usize,
}

impl Drop for ScratchMatrix {
    fn drop(&mut self) {
        let _ = std::fs::remove_file(&self.path);
    }
}

fn create_scratch_matrix() -> anyhow::Result<ScratchMatrix> {
    for _ in 0..1000 {
        let serial = TEMP_COUNTER.fetch_add(1, Ordering::Relaxed);
        let path = std::env::temp_dir().join(format!(
            "phyluce-concatenate-{}-{serial}.matrix",
            std::process::id()
        ));
        match OpenOptions::new()
            .read(true)
            .write(true)
            .create_new(true)
            .open(&path)
        {
            Ok(file) => return Ok(ScratchMatrix { path, file }),
            Err(error) if error.kind() == std::io::ErrorKind::AlreadyExists => continue,
            Err(error) => return Err(error.into()),
        }
    }
    anyhow::bail!("could not allocate a concatenation scratch file")
}

fn collect_metadata(
    files: &[PathBuf],
    input_format: &str,
) -> anyhow::Result<ConcatenationMetadata> {
    let mut taxa = Vec::new();
    let mut taxon_index = HashMap::new();
    let mut charsets = Vec::with_capacity(files.len());
    let mut locus_taxa = Vec::with_capacity(files.len());
    let mut total_length = 0usize;

    for file in files {
        let alignment = load_alignment(file, input_format)
            .with_context(|| format!("loading alignment {}", file.display()))?;
        for row in &alignment.rows {
            if !taxon_index.contains_key(&row.id) {
                taxon_index.insert(row.id.clone(), taxa.len());
                taxa.push(row.id.clone());
            }
        }
        locus_taxa.push(alignment.rows.iter().map(|row| row.id.clone()).collect());
        let stop = total_length
            .checked_add(alignment.nchar())
            .ok_or_else(|| anyhow::anyhow!("concatenated alignment length overflow"))?;
        charsets.push(Charset {
            name: file
                .file_name()
                .and_then(|name| name.to_str())
                .unwrap_or("")
                .to_string(),
            start: total_length,
            stop,
        });
        total_length = stop;
    }
    Ok(ConcatenationMetadata {
        taxa,
        taxon_index,
        charsets,
        locus_taxa,
        total_length,
    })
}

fn initialize_missing_matrix(
    scratch: &mut File,
    taxa: usize,
    total_length: usize,
) -> anyhow::Result<()> {
    let matrix_size = taxa
        .checked_mul(total_length)
        .ok_or_else(|| anyhow::anyhow!("concatenated matrix size overflow"))?;
    let missing = [b'?'; 64 * 1024];
    let mut remaining = matrix_size;
    while remaining > 0 {
        let count = remaining.min(missing.len());
        scratch.write_all(&missing[..count])?;
        remaining -= count;
    }
    scratch.flush()?;
    Ok(())
}

fn populate_matrix(
    scratch: &mut File,
    files: &[PathBuf],
    input_format: &str,
    taxon_index: &HashMap<String, usize>,
    charsets: &[Charset],
    locus_taxa: &[HashSet<String>],
    total_length: usize,
) -> anyhow::Result<()> {
    for ((file, charset), expected_taxa) in files.iter().zip(charsets).zip(locus_taxa) {
        let alignment = load_alignment(file, input_format)
            .with_context(|| format!("loading alignment {}", file.display()))?;
        let expected_length = charset.stop - charset.start;
        anyhow::ensure!(
            alignment.nchar() == expected_length,
            "alignment {} changed between concatenation passes: expected {expected_length} characters, found {}",
            file.display(),
            alignment.nchar()
        );
        let actual_taxa: HashSet<&str> = alignment.rows.iter().map(|row| row.id.as_str()).collect();
        anyhow::ensure!(
            actual_taxa.len() == expected_taxa.len()
                && expected_taxa
                    .iter()
                    .all(|taxon| actual_taxa.contains(taxon.as_str())),
            "alignment {} changed taxon membership between concatenation passes",
            file.display()
        );
        for row in alignment.rows {
            let row_index = *taxon_index.get(&row.id).ok_or_else(|| {
                anyhow::anyhow!(
                    "alignment {} introduced taxon {:?} between concatenation passes",
                    file.display(),
                    row.id
                )
            })?;
            let offset = row_index
                .checked_mul(total_length)
                .and_then(|value| value.checked_add(charset.start))
                .ok_or_else(|| anyhow::anyhow!("concatenated matrix offset overflow"))?;
            scratch.seek(SeekFrom::Start(offset as u64))?;
            scratch.write_all(&row.seq)?;
        }
    }
    scratch.flush()?;
    Ok(())
}

fn copy_taxon_sequence(
    scratch: &mut File,
    output: &mut impl Write,
    taxon_index: usize,
    total_length: usize,
) -> anyhow::Result<()> {
    let offset = taxon_index
        .checked_mul(total_length)
        .ok_or_else(|| anyhow::anyhow!("concatenated matrix offset overflow"))?;
    scratch.seek(SeekFrom::Start(offset as u64))?;
    std::io::copy(&mut scratch.take(total_length as u64), output)?;
    Ok(())
}

fn write_phylip(
    path: &Path,
    scratch: &mut File,
    taxa: &[String],
    total_length: usize,
) -> anyhow::Result<()> {
    let mut output = BufWriter::new(
        File::create(path)
            .with_context(|| format!("creating phylip output file {}", path.display()))?,
    );
    writeln!(output, "{} {}", taxa.len(), total_length)?;
    for (index, taxon) in taxa.iter().enumerate() {
        write!(output, "{} ", safename(taxon))?;
        copy_taxon_sequence(scratch, &mut output, index, total_length)?;
        writeln!(output)?;
    }
    output.flush()?;
    Ok(())
}

fn write_nexus(
    path: &Path,
    scratch: &mut File,
    taxa: &[String],
    total_length: usize,
    charsets: &[Charset],
) -> anyhow::Result<()> {
    let quoted: Vec<String> = taxa.iter().map(|taxon| safename(taxon)).collect();
    let name_length = quoted
        .iter()
        .map(|name| name.chars().count())
        .max()
        .unwrap_or(0);
    let mut output = BufWriter::new(
        File::create(path)
            .with_context(|| format!("creating nexus output file {}", path.display()))?,
    );
    writeln!(output, "#NEXUS\nbegin data;")?;
    writeln!(
        output,
        "dimensions ntax={} nchar={total_length};",
        taxa.len()
    )?;
    writeln!(output, "format datatype=dna missing=? gap=-;\nmatrix")?;
    if total_length > 0 {
        for (index, name) in quoted.iter().enumerate() {
            write!(
                output,
                "{name}{}",
                " ".repeat(name_length + 1 - name.chars().count())
            )?;
            copy_taxon_sequence(scratch, &mut output, index, total_length)?;
            writeln!(output)?;
        }
    }
    write!(output, ";\nend;\n{}", format_sets_block(charsets))?;
    output.flush()?;
    Ok(())
}

pub fn run(
    alignments_dir: &Path,
    input_format: &str,
    output_dir: &Path,
    as_nexus: bool,
    as_phylip: bool,
) -> anyhow::Result<()> {
    anyhow::ensure!(
        as_nexus != as_phylip,
        "exactly one of --nexus or --phylip is required"
    );
    crate::output_path::prepare_output_dir(output_dir)
        .with_context(|| format!("preparing output directory {}", output_dir.display()))?;

    let mut files = find_alignment_files(alignments_dir, input_format)
        .with_context(|| format!("reading alignments directory {}", alignments_dir.display()))?;
    files.sort_by(|left, right| left.file_name().cmp(&right.file_name()));
    let metadata = collect_metadata(&files, input_format)?;
    let mut scratch = create_scratch_matrix()?;
    initialize_missing_matrix(
        &mut scratch.file,
        metadata.taxa.len(),
        metadata.total_length,
    )?;
    populate_matrix(
        &mut scratch.file,
        &files,
        input_format,
        &metadata.taxon_index,
        &metadata.charsets,
        &metadata.locus_taxa,
        metadata.total_length,
    )?;

    let output_name = output_dir
        .file_name()
        .and_then(|name| name.to_str())
        .unwrap_or("concatenated");
    if as_phylip {
        write_phylip(
            &output_dir.join(format!("{output_name}.phylip")),
            &mut scratch.file,
            &metadata.taxa,
            metadata.total_length,
        )?;
        let charsets_path = output_dir.join(format!("{output_name}.charsets"));
        std::fs::write(&charsets_path, format_sets_block(&metadata.charsets))
            .with_context(|| format!("writing charsets file {}", charsets_path.display()))?;
    } else {
        write_nexus(
            &output_dir.join(format!("{output_name}.nexus")),
            &mut scratch.file,
            &metadata.taxa,
            metadata.total_length,
            &metadata.charsets,
        )?;
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn rejects_alignment_changes_between_passes() {
        let root =
            std::env::temp_dir().join(format!("phyluce-concat-consistency-{}", std::process::id()));
        let _ = std::fs::remove_dir_all(&root);
        std::fs::create_dir_all(&root).unwrap();
        let input = root.join("locus.fasta");
        std::fs::write(&input, ">a\nAAAA\n>b\nAAAA\n").unwrap();
        let files = vec![input.clone()];
        let metadata = collect_metadata(&files, "fasta").unwrap();
        let mut scratch = create_scratch_matrix().unwrap();
        initialize_missing_matrix(
            &mut scratch.file,
            metadata.taxa.len(),
            metadata.total_length,
        )
        .unwrap();

        std::fs::write(&input, ">a\nAA\n>b\nAA\n").unwrap();
        let error = populate_matrix(
            &mut scratch.file,
            &files,
            "fasta",
            &metadata.taxon_index,
            &metadata.charsets,
            &metadata.locus_taxa,
            metadata.total_length,
        )
        .unwrap_err();
        assert!(error
            .to_string()
            .contains("changed between concatenation passes"));
        std::fs::remove_dir_all(root).unwrap();
    }
}
