//! CLI wiring for `phyluce align concatenate-alignments`, mirroring
//! `phyluce_align_concatenate_alignments` with a budgeted in-memory or disk-backed matrix.

use std::collections::HashMap;
use std::fs::{File, OpenOptions};
use std::io::{BufReader, BufWriter, Cursor, Read, Seek, SeekFrom, Write};
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
    charsets: Vec<Charset>,
    cached: Option<Vec<PackedLocus>>,
    spool: Option<Spool>,
    total_length: usize,
}

// Internal one-pass snapshot, removed automatically on success or failure.
struct Spool {
    scratch: ScratchMatrix,
    rows: Vec<Vec<DiskFragment>>,
    bytes: u64,
}

impl Spool {
    fn append<'a>(
        &mut self,
        start: usize,
        rows: impl Iterator<Item = (usize, &'a [u8])>,
    ) -> anyhow::Result<()> {
        let mut writer = BufWriter::new(&mut self.scratch.file);
        for (taxon, sequence) in rows {
            if self.rows.len() <= taxon {
                self.rows.resize_with(taxon + 1, Vec::new);
            }
            self.rows[taxon].push(DiskFragment {
                start,
                offset: self.bytes,
                length: sequence.len(),
            });
            writer.write_all(sequence)?;
            self.bytes = self
                .bytes
                .checked_add(sequence.len() as u64)
                .context("snapshot size overflow")?;
        }
        writer.flush()?;
        Ok(())
    }
}

struct DiskFragment {
    start: usize,
    offset: u64,
    length: usize,
}

// Virtual dense matrix over a sparse, append-only snapshot. Missing cells
// never occupy temporary disk space, and no matrix writes are necessary.
struct DiskMatrix {
    spool: Spool,
    width: usize,
    taxa: usize,
    position: u64,
}

impl Read for DiskMatrix {
    fn read(&mut self, bytes: &mut [u8]) -> std::io::Result<usize> {
        let size = self.taxa * self.width;
        let start = usize::try_from(self.position).unwrap_or(usize::MAX);
        let count = bytes.len().min(size.saturating_sub(start));
        bytes[..count].fill(b'?');
        let mut copied = 0;
        while copied < count {
            let offset = start + copied;
            let column = offset % self.width;
            let length = (self.width - column).min(count - copied);
            if let Some(row) = self.spool.rows.get(offset / self.width) {
                let first =
                    row.partition_point(|fragment| fragment.start + fragment.length <= column);
                for fragment in &row[first..] {
                    if fragment.start >= column + length {
                        break;
                    }
                    let left = column.max(fragment.start);
                    let right = (column + length).min(fragment.start + fragment.length);
                    let source = fragment.offset + (left - fragment.start) as u64;
                    let target = &mut bytes[copied + left - column..copied + right - column];
                    #[cfg(unix)]
                    {
                        use std::os::unix::fs::FileExt;
                        self.spool.scratch.file.read_exact_at(target, source)?;
                    }
                    #[cfg(not(unix))]
                    {
                        self.spool.scratch.file.seek(SeekFrom::Start(source))?;
                        self.spool.scratch.file.read_exact(target)?;
                    }
                }
            }
            copied += length;
        }
        self.position += count as u64;
        Ok(count)
    }
}

impl Seek for DiskMatrix {
    fn seek(&mut self, from: SeekFrom) -> std::io::Result<u64> {
        let next = match from {
            SeekFrom::Start(position) => Some(position),
            SeekFrom::Current(offset) => self.position.checked_add_signed(offset),
            SeekFrom::End(offset) => ((self.taxa * self.width) as u64).checked_add_signed(offset),
        }
        .ok_or_else(|| {
            std::io::Error::new(std::io::ErrorKind::InvalidInput, "invalid matrix seek")
        })?;
        self.position = next;
        Ok(next)
    }
}

// One sequence allocation per locus, with global taxon indices in row order.
struct PackedLocus {
    taxa: Box<[usize]>,
    sequence: Box<[u8]>,
    width: usize,
}

struct Fragment {
    start: usize,
    locus: usize,
    offset: usize,
    length: usize,
}

enum Matrix {
    Memory(Cursor<Vec<u8>>),
    Disk(DiskMatrix),
    Cached(CachedMatrix),
}

// Sparse taxon rows retain only observed sequence blocks; missing cells are
// generated during reads. No second, dense matrix allocation is necessary.
struct CachedMatrix {
    rows: Vec<Vec<Fragment>>,
    loci: Vec<PackedLocus>,
    width: usize,
    position: u64,
}

impl CachedMatrix {
    fn new(metadata: &ConcatenationMetadata, loci: Vec<PackedLocus>) -> Self {
        let mut counts = vec![0; metadata.taxa.len()];
        for locus in &loci {
            for &taxon in &locus.taxa {
                counts[taxon] += 1;
            }
        }
        let mut rows: Vec<Vec<Fragment>> = counts.into_iter().map(Vec::with_capacity).collect();
        for (index, (locus, charset)) in loci.iter().zip(&metadata.charsets).enumerate() {
            for (row, &taxon) in locus.taxa.iter().enumerate() {
                rows[taxon].push(Fragment {
                    start: charset.start,
                    locus: index,
                    offset: row * locus.width,
                    length: locus.width,
                });
            }
        }
        Self {
            rows,
            loci,
            width: metadata.total_length,
            position: 0,
        }
    }
}

impl Read for CachedMatrix {
    fn read(&mut self, bytes: &mut [u8]) -> std::io::Result<usize> {
        let size = self.rows.len() * self.width;
        let start = usize::try_from(self.position).unwrap_or(usize::MAX);
        let count = bytes.len().min(size.saturating_sub(start));
        bytes[..count].fill(b'?');
        let mut copied = 0;
        while copied < count {
            let offset = start + copied;
            let column = offset % self.width;
            let length = (self.width - column).min(count - copied);
            let row = &self.rows[offset / self.width];
            let first = row.partition_point(|fragment| fragment.start + fragment.length <= column);
            for fragment in &row[first..] {
                let begin = fragment.start;
                if begin >= column + length {
                    break;
                }
                let left = column.max(begin);
                let right = (column + length).min(begin + fragment.length);
                let sequence = &self.loci[fragment.locus].sequence;
                bytes[copied + left - column..copied + right - column].copy_from_slice(
                    &sequence[fragment.offset + left - begin..fragment.offset + right - begin],
                );
            }
            copied += length;
        }
        self.position += count as u64;
        Ok(count)
    }
}

impl Seek for CachedMatrix {
    fn seek(&mut self, from: SeekFrom) -> std::io::Result<u64> {
        let next = match from {
            SeekFrom::Start(position) => Some(position),
            SeekFrom::Current(offset) => self.position.checked_add_signed(offset),
            SeekFrom::End(offset) => {
                ((self.rows.len() * self.width) as u64).checked_add_signed(offset)
            }
        }
        .ok_or_else(|| {
            std::io::Error::new(std::io::ErrorKind::InvalidInput, "invalid matrix seek")
        })?;
        self.position = next;
        Ok(next)
    }
}

impl Read for Matrix {
    fn read(&mut self, bytes: &mut [u8]) -> std::io::Result<usize> {
        match self {
            Self::Memory(matrix) => matrix.read(bytes),
            Self::Cached(matrix) => matrix.read(bytes),
            Self::Disk(matrix) => matrix.read(bytes),
        }
    }
}
impl Seek for Matrix {
    fn seek(&mut self, position: SeekFrom) -> std::io::Result<u64> {
        match self {
            Self::Memory(matrix) => matrix.seek(position),
            Self::Cached(matrix) => matrix.seek(position),
            Self::Disk(matrix) => matrix.seek(position),
        }
    }
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
    budget: usize,
) -> anyhow::Result<ConcatenationMetadata> {
    let mut taxa = Vec::new();
    let mut taxon_index = HashMap::new();
    let mut charsets = Vec::with_capacity(files.len());
    let mut spool: Option<Spool> = None;
    let mut total_length = 0usize;
    let mut cached = Some(Vec::new());
    let mut cached_size = 0usize;

    for file in files {
        let alignment = load_alignment(file, input_format)
            .with_context(|| format!("loading alignment {}", file.display()))?;
        for row in &alignment.rows {
            if !taxon_index.contains_key(&row.id) {
                taxon_index.insert(row.id.clone(), taxa.len());
                taxa.push(row.id.clone());
            }
        }
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
        let _matrix_size = taxa
            .len()
            .checked_mul(total_length)
            .context("concatenated matrix size overflow")?;
        if let Some(alignments) = &mut cached {
            let width = alignment.nchar();
            let count = alignment.rows.len();
            // Include the packed data, taxon indices and eventual row fragments.
            let bytes = count.checked_mul(width);
            let next = count
                .checked_mul(std::mem::size_of::<usize>() + std::mem::size_of::<Fragment>())
                .and_then(|n| n.checked_add(bytes?))
                .and_then(|n| n.checked_add(std::mem::size_of::<PackedLocus>()))
                .and_then(|n| cached_size.checked_add(n));
            if next.is_some_and(|size| size <= budget) {
                let mut sequence = Vec::with_capacity(bytes.unwrap());
                let mut indices = Vec::with_capacity(count);
                for row in &alignment.rows {
                    indices.push(taxon_index[&row.id]);
                    sequence.extend_from_slice(&row.seq);
                }
                cached_size = next.unwrap();
                alignments.push(PackedLocus {
                    taxa: indices.into_boxed_slice(),
                    sequence: sequence.into_boxed_slice(),
                    width,
                });
            } else {
                let mut disk = Spool {
                    scratch: create_scratch_matrix()?,
                    rows: Vec::new(),
                    bytes: 0,
                };
                for (locus, charset) in cached.take().unwrap().into_iter().zip(&charsets) {
                    disk.append(
                        charset.start,
                        locus.taxa.iter().enumerate().map(|(row, &taxon)| {
                            (
                                taxon,
                                &locus.sequence[row * locus.width..(row + 1) * locus.width],
                            )
                        }),
                    )?;
                }
                spool = Some(disk);
            }
        }
        if let Some(disk) = &mut spool {
            disk.append(
                charsets.last().unwrap().start,
                alignment
                    .rows
                    .iter()
                    .map(|row| (taxon_index[&row.id], row.seq.as_slice())),
            )?;
        }
    }
    Ok(ConcatenationMetadata {
        taxa,
        charsets,
        cached,
        spool,
        total_length,
    })
}

fn spilled_matrix(metadata: &mut ConcatenationMetadata, budget: usize) -> anyhow::Result<Matrix> {
    let size = metadata.taxa.len() * metadata.total_length;
    let mut bytes = Vec::new();
    if budget > 0 && size <= budget && bytes.try_reserve_exact(size).is_ok() {
        bytes.resize(size, b'?');
        let mut spool = metadata
            .spool
            .take()
            .context("missing concatenation snapshot")?;
        // Revisit the append order for sequential disk reads. These references
        // are temporary and do not duplicate sequence data.
        let mut fragments: Vec<_> = spool
            .rows
            .iter()
            .enumerate()
            .flat_map(|(taxon, row)| row.iter().map(move |fragment| (taxon, fragment)))
            .collect();
        fragments.sort_unstable_by_key(|(_, fragment)| fragment.offset);
        spool.scratch.file.seek(SeekFrom::Start(0))?;
        let mut reader = BufReader::new(&mut spool.scratch.file);
        for (taxon, fragment) in fragments {
            let offset = taxon * metadata.total_length + fragment.start;
            reader.read_exact(&mut bytes[offset..offset + fragment.length])?;
        }
        Ok(Matrix::Memory(Cursor::new(bytes)))
    } else {
        disk_matrix(metadata)
    }
}

fn disk_matrix(metadata: &mut ConcatenationMetadata) -> anyhow::Result<Matrix> {
    Ok(Matrix::Disk(DiskMatrix {
        spool: metadata
            .spool
            .take()
            .context("missing concatenation snapshot")?,
        width: metadata.total_length,
        taxa: metadata.taxa.len(),
        position: 0,
    }))
}

fn copy_taxon_sequence(
    scratch: &mut (impl Read + Seek),
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
    scratch: &mut (impl Read + Seek),
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
    scratch: &mut (impl Read + Seek),
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
    memory_limit_mb: usize,
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
    let budget = memory_limit_mb
        .checked_mul(1024 * 1024)
        .context("memory budget overflow")?;
    let mut metadata = collect_metadata(&files, input_format, budget)?;
    let size = metadata
        .taxa
        .len()
        .checked_mul(metadata.total_length)
        .context("concatenated matrix size overflow")?;
    let mut matrix = if let Some(cached) = metadata.cached.take() {
        Matrix::Cached(CachedMatrix::new(&metadata, cached))
    } else {
        spilled_matrix(&mut metadata, budget)?
    };
    tracing::info!(
        matrix_bytes = size,
        cached_inputs = matches!(&matrix, Matrix::Cached(_)),
        in_memory = !matches!(&matrix, Matrix::Disk(_)),
        "concatenation storage"
    );

    let output_name = output_dir
        .file_name()
        .and_then(|name| name.to_str())
        .unwrap_or("concatenated");
    if as_phylip {
        write_phylip(
            &output_dir.join(format!("{output_name}.phylip")),
            &mut matrix,
            &metadata.taxa,
            metadata.total_length,
        )?;
        let charsets_path = output_dir.join(format!("{output_name}.charsets"));
        std::fs::write(&charsets_path, format_sets_block(&metadata.charsets))
            .with_context(|| format!("writing charsets file {}", charsets_path.display()))?;
    } else {
        write_nexus(
            &output_dir.join(format!("{output_name}.nexus")),
            &mut matrix,
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
    fn cached_and_disk_matrices_match_with_missing_taxa() {
        let root = std::env::temp_dir().join(format!("phyluce-cache-test-{}", std::process::id()));
        std::fs::create_dir_all(&root).unwrap();
        let files = vec![root.join("a.fasta"), root.join("b.fasta")];
        std::fs::write(&files[0], ">z\nAC\n>a\nGT\n").unwrap();
        std::fs::write(&files[1], ">b\nAAA\n>z\nCCC\n").unwrap();
        let mut metadata = collect_metadata(&files, "fasta", 1024 * 1024).unwrap();
        let cached = metadata.cached.take().unwrap();
        let mut sparse = CachedMatrix::new(&metadata, cached);
        let expected = b"ACCCCGT?????AAA";
        for start in 0..=expected.len() {
            sparse.seek(SeekFrom::Start(start as u64)).unwrap();
            let mut bytes = Vec::new();
            let mut buffer = [0; 3];
            loop {
                let n = sparse.read(&mut buffer).unwrap();
                if n == 0 {
                    break;
                }
                bytes.extend_from_slice(&buffer[..n]);
            }
            assert_eq!(bytes, expected[start..]);
        }
        assert!(sparse.seek(SeekFrom::End(-16)).is_err());
        sparse.seek(SeekFrom::End(10)).unwrap();
        assert_eq!(sparse.read(&mut [0; 1]).unwrap(), 0);
        for budget in [0, 15, 150, 1024 * 1024] {
            let mut metadata = collect_metadata(&files, "fasta", budget).unwrap();
            assert_eq!(metadata.cached.is_some(), budget > 150);
            let cached = metadata.cached.take();
            if cached.is_some() {
                for file in &files {
                    std::fs::remove_file(file).unwrap();
                }
            }
            let mut matrix = if let Some(cached) = cached {
                Matrix::Cached(CachedMatrix::new(&metadata, cached))
            } else {
                spilled_matrix(&mut metadata, budget).unwrap()
            };
            matrix.seek(SeekFrom::Start(0)).unwrap();
            let mut data = Vec::new();
            matrix.read_to_end(&mut data).unwrap();
            assert_eq!(data, b"ACCCCGT?????AAA");
            for start in 0..=data.len() {
                matrix.seek(SeekFrom::Start(start as u64)).unwrap();
                let mut small = [0; 3];
                let count = matrix.read(&mut small).unwrap();
                assert_eq!(&small[..count], &data[start..(start + 3).min(data.len())]);
            }
            assert!(matrix.seek(SeekFrom::End(-16)).is_err());
        }
        std::fs::remove_dir_all(root).unwrap();
    }

    #[test]
    fn snapshot_survives_input_removal_and_rejects_truncation() {
        let root =
            std::env::temp_dir().join(format!("phyluce-concat-consistency-{}", std::process::id()));
        let _ = std::fs::remove_dir_all(&root);
        std::fs::create_dir_all(&root).unwrap();
        let input = root.join("locus.fasta");
        std::fs::write(&input, ">a\nAAAA\n>b\nAAAA\n").unwrap();
        let files = vec![input.clone()];
        let mut metadata = collect_metadata(&files, "fasta", 0).unwrap();
        let mut broken = collect_metadata(&files, "fasta", 0).unwrap();
        let snapshot_path = metadata.spool.as_ref().unwrap().scratch.path.clone();
        std::fs::remove_file(input).unwrap();
        let mut matrix = disk_matrix(&mut metadata).unwrap();
        let mut bytes = Vec::new();
        matrix.read_to_end(&mut bytes).unwrap();
        assert_eq!(bytes, b"AAAAAAAA");
        assert_eq!(std::fs::metadata(&snapshot_path).unwrap().len(), 8);
        drop(matrix);
        assert!(!snapshot_path.exists());
        broken
            .spool
            .as_mut()
            .unwrap()
            .scratch
            .file
            .set_len(2)
            .unwrap();
        let broken_path = broken.spool.as_ref().unwrap().scratch.path.clone();
        let mut matrix = disk_matrix(&mut broken).unwrap();
        assert_eq!(
            matrix.read_to_end(&mut Vec::new()).unwrap_err().kind(),
            std::io::ErrorKind::UnexpectedEof
        );
        drop(matrix);
        assert!(!broken_path.exists());
        std::fs::remove_dir_all(root).unwrap();
    }
}
