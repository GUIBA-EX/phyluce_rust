//! phyluce-io: FASTA/FASTQ readers and writers shared across phyluce commands.
//!
//! This first slice only covers FASTA reading (plain + gzip), enough to back
//! `phyluce io validate-fasta` and `phyluce assembly get-fasta-lengths`.

use std::fs::File;
use std::io::{self, BufRead, BufReader, Read};
use std::path::Path;

#[derive(Debug, thiserror::Error)]
pub enum FastaError {
    #[error("{0}")]
    Io(#[from] io::Error),
    #[error("line {line}: sequence data before any header")]
    DataBeforeHeader { line: usize },
    #[error("empty file")]
    Empty,
}

/// One FASTA record: header split into `id` (first whitespace-delimited
/// token, `>` stripped) and the remainder of the header line as `description`.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct FastaRecord {
    pub id: String,
    pub description: String,
    pub sequence: String,
}

pub mod fastq;
pub mod lastz;
pub mod sql;
pub mod twobit;

pub(crate) fn open_maybe_gz(path: &Path) -> io::Result<Box<dyn BufRead>> {
    let f = File::open(path)?;
    if path.extension().and_then(|e| e.to_str()) == Some("gz") {
        Ok(Box::new(BufReader::new(flate2::read::GzDecoder::new(f))))
    } else {
        Ok(Box::new(BufReader::new(f)))
    }
}

/// Parse a full FASTA file into records. Sequence lines are concatenated
/// verbatim (no whitespace stripped beyond line endings) to make roundtrips
/// exact; use `.sequence.len()` for length-only use cases.
pub fn read_fasta(path: &Path) -> Result<Vec<FastaRecord>, FastaError> {
    let reader = open_maybe_gz(path)?;
    read_fasta_reader(reader)
}

/// Parse FASTA records from an existing buffered reader.
pub fn read_fasta_reader(mut reader: impl BufRead) -> Result<Vec<FastaRecord>, FastaError> {
    let mut records = Vec::new();
    let mut current: Option<(String, String, String)> = None;
    let mut line = String::new();
    let mut line_number = 0usize;

    loop {
        line.clear();
        if reader.read_line(&mut line)? == 0 {
            break;
        }
        line_number += 1;
        // One `trim_end_matches` instead of two `ends_with`+`pop` calls
        // followed by a *second*, full `.trim()` re-scan of the sequence
        // branch below (once to check emptiness, once again to get the
        // trimmed slice to push): trailing `\n`/`\r` only ever needs
        // stripping once, up front.
        let trimmed = line.trim_end_matches(['\n', '\r']);
        if let Some(header) = trimmed.strip_prefix('>') {
            if let Some((id, description, sequence)) = current.take() {
                records.push(FastaRecord {
                    id,
                    description,
                    sequence,
                });
            }
            let id = header.split_whitespace().next().unwrap_or("").to_string();
            current = Some((id, header.to_string(), String::new()));
        } else {
            let seq_line = trimmed.trim();
            if !seq_line.is_empty() {
                match &mut current {
                    Some((_, _, seq)) => seq.push_str(seq_line),
                    None => return Err(FastaError::DataBeforeHeader { line: line_number }),
                }
            }
        }
    }
    if let Some((id, description, sequence)) = current {
        records.push(FastaRecord {
            id,
            description,
            sequence,
        });
    } else if records.is_empty() {
        return Err(FastaError::Empty);
    }
    Ok(records)
}

/// Per-record sequence lengths, mirroring the legacy `fasta_iter` generator
/// used by `phyluce_assembly_get_fasta_lengths` (sum of stripped sequence
/// line lengths following each header line).
pub fn fasta_lengths(path: &Path) -> Result<Vec<usize>, FastaError> {
    Ok(read_fasta(path)?
        .into_iter()
        .map(|r| r.sequence.len())
        .collect())
}

/// Validation issue found while checking a FASTA file's structure.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct ValidationIssue {
    pub line: usize,
    pub message: String,
}

/// Validate that a file is well-formed FASTA: starts with `>`, every record
/// has a non-empty id, non-empty sequence, and only IUPAC nucleotide/amino
/// acid or gap/ambiguity characters.
pub fn validate_fasta(path: &Path) -> Result<Vec<ValidationIssue>, FastaError> {
    let reader = open_maybe_gz(path)?;
    let mut issues = Vec::new();
    let mut seen_header = false;
    let mut current_id: Option<String> = None;
    let mut current_len = 0usize;
    let mut current_header_line = 0usize;

    const VALID_CHARS: &str = "ACGTUNRYSWKMBDHVacgtunryswkmbdhv-.*XZBJxzbj";

    for (i, line) in reader.lines().enumerate() {
        let line = line?;
        let lineno = i + 1;
        if let Some(header) = line.strip_prefix('>') {
            if seen_header && current_len == 0 {
                issues.push(ValidationIssue {
                    line: current_header_line,
                    message: format!(
                        "record '{}' has an empty sequence",
                        current_id.clone().unwrap_or_default()
                    ),
                });
            }
            seen_header = true;
            current_header_line = lineno;
            current_len = 0;
            let id = header.split_whitespace().next().unwrap_or("");
            if id.is_empty() {
                issues.push(ValidationIssue {
                    line: lineno,
                    message: "header has no identifier".to_string(),
                });
            }
            current_id = Some(id.to_string());
        } else if !line.trim().is_empty() {
            if !seen_header {
                issues.push(ValidationIssue {
                    line: lineno,
                    message: "sequence data before any header".to_string(),
                });
                continue;
            }
            for c in line.trim().chars() {
                if !VALID_CHARS.contains(c) {
                    issues.push(ValidationIssue {
                        line: lineno,
                        message: format!("invalid sequence character '{}'", c),
                    });
                    break;
                }
            }
            current_len += line.trim().len();
        }
    }
    if !seen_header {
        issues.push(ValidationIssue {
            line: 0,
            message: "file contains no FASTA headers".to_string(),
        });
    } else if current_len == 0 {
        issues.push(ValidationIssue {
            line: current_header_line,
            message: format!(
                "record '{}' has an empty sequence",
                current_id.unwrap_or_default()
            ),
        });
    }
    Ok(issues)
}

/// Read raw bytes from a maybe-gzipped file (used by callers that want to
/// hash/compare content wholesale rather than parse it).
pub fn read_all(path: &Path) -> io::Result<Vec<u8>> {
    let mut buf = Vec::new();
    open_maybe_gz(path)?.read_to_end(&mut buf)?;
    Ok(buf)
}

/// Write one FASTA record wrapped at 60 characters per line, matching
/// Biopython's `SeqRecord.format("fasta")` (`Bio.SeqIO.FastaIO`'s default
/// wrap width). `id` is written verbatim as the whole header line (no
/// separate description is appended) since several phyluce commands stuff
/// the full `id |extra` header text into `record.id` itself.
pub fn write_fasta_record<W: std::io::Write>(
    writer: &mut W,
    id: &str,
    sequence: &str,
) -> io::Result<()> {
    writeln!(writer, ">{id}")?;
    for chunk in sequence.as_bytes().chunks(60) {
        writer.write_all(chunk)?;
        writer.write_all(b"\n")?;
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::Write;

    fn write_temp(name: &str, contents: &str) -> std::path::PathBuf {
        let dir = std::env::temp_dir().join("phyluce-io-tests");
        std::fs::create_dir_all(&dir).unwrap();
        let path = dir.join(name);
        let mut f = File::create(&path).unwrap();
        f.write_all(contents.as_bytes()).unwrap();
        path
    }

    // Ad hoc benchmark for `read_fasta`'s per-line String allocation.
    // Run with:
    //   cargo +stable test --release -p phyluce-io --lib -- --ignored --nocapture bench_
    //
    // Checked against rust-bio's `bio::io::fasta`/`fastq` readers (a
    // temporary dev-dependency, removed after measuring): on this same
    // benchmark, rust-bio parsed ~2x faster for FASTA and ~1.8x faster for
    // FASTQ (both still >1 GB/s here either way). Not adopted: `bio` pulls
    // in 83 transitive crates (ndarray, nalgebra, statrs, petgraph, rand,
    // ...) for functionality this project only needs a sliver of, and
    // phyluce-io is a foundational dependency of every command in this
    // workspace -- that cost would be paid everywhere for a component
    // that real UCE pipelines don't spend much time in (LASTZ/MAFFT/
    // SPAdes/etc. subprocess time dominates, not our own file parsing).
    //
    // Hand-optimized this reader afterward instead (no new dependency):
    // the old version re-scanned each sequence line up to 3x (`ends_with`
    // + `pop` for `\n`, `\r`, then `.trim().is_empty()`, then `.trim()`
    // again for the `push_str` argument) where rust-bio's equivalent loop
    // does one `trim_end()`. Collapsing that to a single
    // `trim_end_matches` + one `.trim()` closed part of the gap (~1.1-1.2
    // GB/s -> ~1.35-1.4 GB/s here), and the same fix applied to
    // `fastq::fastq_lengths`/`fastq_record_count` (which used to allocate
    // and UTF-8-validate a fresh `String` per line via `BufRead::lines()`
    // for every line, including the 3-in-4 they never read) got
    // `fastq_lengths` from ~1.65 GB/s to ~1.95 GB/s and
    // `fastq_record_count` -- which never even looks at line content -- to
    // ~5.5 GB/s. Still short of rust-bio's ~2.4/~3.1 GB/s on the same
    // data; closing the rest would mean working in raw bytes throughout
    // instead of `String` (a bigger change to `FastaRecord`'s public
    // shape), which isn't justified without evidence this is an actual
    // bottleneck in a real run.
    #[test]
    #[ignore]
    fn bench_read_fasta_large_file() {
        // ~5000 records * ~500bp, wrapped at 70 cols -- resembles a
        // mid-size contig assembly FASTA.
        let dir = std::env::temp_dir().join("phyluce-io-bench");
        std::fs::create_dir_all(&dir).unwrap();
        let path = dir.join("bench.fasta");
        {
            let mut f = std::io::BufWriter::new(File::create(&path).unwrap());
            let bases = b"ACGT";
            for i in 0..5_000usize {
                writeln!(f, ">NODE_{i}_length_500_cov_10.3").unwrap();
                for line in 0..8 {
                    let mut buf = [0u8; 70];
                    for (j, b) in buf.iter_mut().enumerate() {
                        *b = bases[(i + line + j) % 4];
                    }
                    f.write_all(&buf).unwrap();
                    f.write_all(b"\n").unwrap();
                }
            }
        }
        let file_len = std::fs::metadata(&path).unwrap().len();

        let start = std::time::Instant::now();
        let records = read_fasta(&path).unwrap();
        let elapsed = start.elapsed();

        eprintln!(
            "[bench] read_fasta: {} bytes / {} records in {:?} ({:.1} MB/s)",
            file_len,
            records.len(),
            elapsed,
            (file_len as f64 / 1_000_000.0) / elapsed.as_secs_f64()
        );
        std::fs::remove_file(&path).ok();
    }

    #[test]
    fn writes_fasta_wrapped_at_60() {
        let seq = "A".repeat(65);
        let mut out = Vec::new();
        write_fasta_record(&mut out, "uce-1", &seq).unwrap();
        let text = String::from_utf8(out).unwrap();
        let expected = format!(">uce-1\n{}\n{}\n", "A".repeat(60), "A".repeat(5));
        assert_eq!(text, expected);
    }

    #[test]
    fn reads_basic_fasta() {
        let path = write_temp(
            "basic.fasta",
            ">uce-1 |taxon\nACGT\nACGT\n>uce-2 |taxon\nAC\n",
        );
        let records = read_fasta(&path).unwrap();
        assert_eq!(records.len(), 2);
        assert_eq!(records[0].id, "uce-1");
        assert_eq!(records[0].sequence, "ACGTACGT");
        assert_eq!(records[1].sequence, "AC");
    }

    #[test]
    fn lengths_match_record_sequence_len() {
        let path = write_temp("lengths.fasta", ">a\nACGT\n>b\nACGTACGT\n");
        let lengths = fasta_lengths(&path).unwrap();
        assert_eq!(lengths, vec![4, 8]);
    }

    #[test]
    fn validate_flags_empty_sequence() {
        let path = write_temp("empty_seq.fasta", ">a\n>b\nACGT\n");
        let issues = validate_fasta(&path).unwrap();
        assert!(issues.iter().any(|i| i.message.contains("empty sequence")));
    }

    #[test]
    fn validate_flags_bad_characters() {
        let path = write_temp("bad_chars.fasta", ">a\nACGT123\n");
        let issues = validate_fasta(&path).unwrap();
        assert!(issues
            .iter()
            .any(|i| i.message.contains("invalid sequence character")));
    }

    #[test]
    fn validate_clean_file_has_no_issues() {
        let path = write_temp("clean.fasta", ">a\nACGT\n>b\nACGT\n");
        let issues = validate_fasta(&path).unwrap();
        assert!(issues.is_empty());
    }
}
