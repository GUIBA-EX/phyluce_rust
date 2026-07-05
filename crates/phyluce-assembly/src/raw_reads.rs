//! Raw FASTQ/FASTA read discovery, mirroring `phyluce/raw_reads.py`: find
//! R1/R2/singleton files in a sample directory and classify them.

use std::path::{Path, PathBuf};

#[derive(Debug, thiserror::Error)]
pub enum RawReadsError {
    #[error("{0}")]
    Io(#[from] std::io::Error),
    #[error("there are both fasta and fastq files in {0}")]
    MixedFastaFastq(String),
    #[error("there are no appropriate files in {0}")]
    NoFiles(String),
    #[error("there appear to be multiple files for R1/R2/Singleton reads")]
    MultipleFilesForRead,
    #[error("files are of different types (e.g. gzip and fastq)")]
    MixedExtensions,
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub enum FileKind {
    Fastq,
    Fasta,
}

#[derive(Debug, Clone, Default)]
pub struct ReadSet {
    pub kind: Option<FileKind>,
    pub r1: Option<PathBuf>,
    pub r2: Option<PathBuf>,
    pub singleton: Option<PathBuf>,
    pub gzip: bool,
}

const FASTQ_EXTENSIONS: &[&str] = &[
    ".fastq.gz",
    ".fastq.gzip",
    ".fq.gz",
    ".fq.gzip",
    ".fq",
    ".fastq",
];
const FASTA_EXTENSIONS: &[&str] = &[
    ".fasta.gz",
    ".fasta.gzip",
    ".fa.gz",
    ".fa.gzip",
    ".fa",
    ".fasta",
];

fn find_by_extensions(
    dir: &Path,
    subfolder: &str,
    extensions: &[&str],
) -> std::io::Result<Vec<PathBuf>> {
    let search_dir = if subfolder.is_empty() {
        dir.to_path_buf()
    } else {
        dir.join(subfolder)
    };
    let mut out = Vec::new();
    if !search_dir.is_dir() {
        return Ok(out);
    }
    for entry in std::fs::read_dir(&search_dir)? {
        let path = entry?.path();
        if let Some(name) = path.file_name().and_then(|n| n.to_str()) {
            if extensions.iter().any(|ext| name.ends_with(ext)) {
                out.push(path);
            }
        }
    }
    out.sort();
    Ok(out)
}

/// Classify a filename as R1 / R2 / singleton-or-unpaired using the same
/// regex family as the rest of phyluce:
/// `(?:.*)[_-](?:READ|Read|R)(\d)*[_-]*(singleton|unpaired)*(?:.*)`.
/// Returns `None` if the filename doesn't match any recognized pattern.
fn classify_read(fname: &str) -> Option<&'static str> {
    // A small hand-rolled matcher instead of pulling in `regex` for this
    // one alternation: scan for a `[_-](READ|Read|R)` marker, then check
    // what immediately follows it.
    let bytes = fname.as_bytes();
    for i in 0..bytes.len() {
        if bytes[i] != b'_' && bytes[i] != b'-' {
            continue;
        }
        let rest = &fname[i + 1..];
        for marker in ["READ", "Read", "R"] {
            if let Some(after) = rest.strip_prefix(marker) {
                let digits: String = after.chars().take_while(|c| c.is_ascii_digit()).collect();
                let after_digits = &after[digits.len()..];
                if digits == "1" {
                    return Some("r1");
                } else if digits == "2" {
                    return Some("r2");
                }
                let trimmed = after_digits.trim_start_matches(['_', '-']);
                if trimmed.starts_with("singleton") || trimmed.starts_with("unpaired") {
                    return Some("singleton");
                }
            }
        }
    }
    None
}

/// Mirrors `get_input_files`: locate and classify every FASTQ/FASTA file
/// in `dir`/`subfolder`.
pub fn get_input_files(dir: &Path, subfolder: &str) -> Result<ReadSet, RawReadsError> {
    let fastq_files = find_by_extensions(dir, subfolder, FASTQ_EXTENSIONS)?;
    let fasta_files = find_by_extensions(dir, subfolder, FASTA_EXTENSIONS)?;

    if !fastq_files.is_empty() && !fasta_files.is_empty() {
        return Err(RawReadsError::MixedFastaFastq(dir.display().to_string()));
    }
    let (files, kind) = if !fastq_files.is_empty() {
        (fastq_files, FileKind::Fastq)
    } else if !fasta_files.is_empty() {
        (fasta_files, FileKind::Fasta)
    } else {
        return Err(RawReadsError::NoFiles(dir.display().to_string()));
    };

    let mut set = ReadSet {
        kind: Some(kind),
        ..Default::default()
    };
    let mut extensions = std::collections::HashSet::new();
    for f in &files {
        let fname = f.file_name().and_then(|n| n.to_str()).unwrap_or("");
        if let Some(dot) = fname.rfind('.') {
            extensions.insert(fname[dot..].to_string());
        }
        match classify_read(fname) {
            Some("r1") => {
                if set.r1.is_some() {
                    return Err(RawReadsError::MultipleFilesForRead);
                }
                set.r1 = Some(f.clone());
            }
            Some("r2") => {
                if set.r2.is_some() {
                    return Err(RawReadsError::MultipleFilesForRead);
                }
                set.r2 = Some(f.clone());
            }
            Some("singleton") => {
                if set.singleton.is_some() {
                    return Err(RawReadsError::MultipleFilesForRead);
                }
                set.singleton = Some(f.clone());
            }
            _ => {}
        }
    }
    if extensions.len() != 1 {
        return Err(RawReadsError::MixedExtensions);
    }
    let ext = extensions.iter().next().unwrap();
    set.gzip = ext.ends_with(".gz") || ext.ends_with(".gzip");
    Ok(set)
}

#[cfg(test)]
mod tests {
    use super::*;

    fn touch(dir: &Path, name: &str) {
        std::fs::write(dir.join(name), b"x").unwrap();
    }

    #[test]
    fn classifies_r1_r2_singleton() {
        assert_eq!(classify_read("sample-READ1.fastq.gz"), Some("r1"));
        assert_eq!(classify_read("sample-READ2.fastq.gz"), Some("r2"));
        assert_eq!(
            classify_read("sample-READ-singleton.fastq.gz"),
            Some("singleton")
        );
        assert_eq!(classify_read("sample_R1_001.fastq.gz"), Some("r1"));
    }

    #[test]
    fn discovers_and_classifies_fastq_set() {
        let dir = std::env::temp_dir().join(format!("phyluce-raw-reads-{}", std::process::id()));
        std::fs::create_dir_all(&dir).unwrap();
        touch(&dir, "sample-READ1.fastq.gz");
        touch(&dir, "sample-READ2.fastq.gz");
        touch(&dir, "sample-READ-singleton.fastq.gz");

        let set = get_input_files(&dir, "").unwrap();
        assert_eq!(set.kind, Some(FileKind::Fastq));
        assert!(set.r1.is_some());
        assert!(set.r2.is_some());
        assert!(set.singleton.is_some());
        assert!(set.gzip);
    }

    #[test]
    fn rejects_mixed_fasta_and_fastq() {
        let dir =
            std::env::temp_dir().join(format!("phyluce-raw-reads-mixed-{}", std::process::id()));
        std::fs::create_dir_all(&dir).unwrap();
        touch(&dir, "sample-READ1.fastq.gz");
        touch(&dir, "sample-READ1.fasta");
        let err = get_input_files(&dir, "").unwrap_err();
        assert!(matches!(err, RawReadsError::MixedFastaFastq(_)));
    }
}
