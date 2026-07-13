//! MAFFT wrapper mirroring `phyluce/mafft.py`'s `Align.run_alignment`.

use std::io::Write as _;
use std::path::{Path, PathBuf};
use std::sync::atomic::{AtomicU64, Ordering};

use phyluce_external::ExternalCommand;
use phyluce_io::read_fasta_reader;

use crate::{Alignment, AlignmentRow};

static TEMP_COUNTER: AtomicU64 = AtomicU64::new(0);

struct UniqueTempFile {
    file: std::fs::File,
    path: PathBuf,
}

impl UniqueTempFile {
    fn new() -> std::io::Result<Self> {
        loop {
            let serial = TEMP_COUNTER.fetch_add(1, Ordering::Relaxed);
            let path = std::env::temp_dir().join(format!(
                "phyluce-mafft-{}-{serial}.fasta",
                std::process::id()
            ));
            match std::fs::OpenOptions::new()
                .write(true)
                .create_new(true)
                .open(&path)
            {
                Ok(file) => return Ok(Self { file, path }),
                Err(error) if error.kind() == std::io::ErrorKind::AlreadyExists => continue,
                Err(error) => return Err(error),
            }
        }
    }

    fn path(&self) -> &Path {
        &self.path
    }
}

impl std::io::Write for UniqueTempFile {
    fn write(&mut self, buffer: &[u8]) -> std::io::Result<usize> {
        self.file.write(buffer)
    }

    fn flush(&mut self) -> std::io::Result<()> {
        self.file.flush()
    }
}

impl Drop for UniqueTempFile {
    fn drop(&mut self) {
        let _ = std::fs::remove_file(&self.path);
    }
}

#[derive(Debug, thiserror::Error)]
pub enum MafftError {
    #[error("{0}")]
    Io(#[from] std::io::Error),
    #[error("{0}")]
    External(#[from] phyluce_external::ExternalError),
    #[error("{0}")]
    Fasta(#[from] phyluce_io::FastaError),
    #[error("{0}")]
    Alignment(#[from] crate::AlignmentError),
}

/// Run MAFFT (`--adjustdirection --maxiterate 1000`) on the given
/// (id, sequence) records and return the resulting alignment. Mirrors
/// `Align.run_alignment`: MAFFT may prepend `_R_` to the id of any
/// sequence it reverse-complemented to fit the alignment -- this is
/// intentionally NOT stripped here (the legacy code doesn't either;
/// several downstream commands strip it explicitly when they need to).
pub fn run_mafft(mafft_bin: &str, records: &[(String, String)]) -> Result<Alignment, MafftError> {
    let mut input = UniqueTempFile::new()?;
    for (id, seq) in records {
        writeln!(input, ">{id}")?;
        writeln!(input, "{seq}")?;
    }
    input.flush()?;

    let report = ExternalCommand::new(mafft_bin)
        .args([
            "--adjustdirection".to_string(),
            "--maxiterate".to_string(),
            "1000".to_string(),
            input.path().to_string_lossy().to_string(),
        ])
        .run()?;
    let parsed = read_fasta_reader(std::io::Cursor::new(report.stdout.as_bytes()))?;

    let alignment = Alignment {
        rows: parsed
            .into_iter()
            .map(|r| AlignmentRow {
                id: r.id,
                seq: r.sequence.into_bytes(),
            })
            .collect(),
    };
    alignment.validate()?;
    Ok(alignment)
}

pub fn mafft_binary_path(
    cfg: &phyluce_config::PhyluceConfig,
) -> Result<String, phyluce_config::ConfigError> {
    cfg.get_user_path("binaries", "mafft")
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn temporary_inputs_are_unique_and_removed_on_drop() {
        let first = UniqueTempFile::new().unwrap();
        let first_path = first.path().to_path_buf();
        let second = UniqueTempFile::new().unwrap();
        let second_path = second.path().to_path_buf();
        assert_ne!(first_path, second_path);
        assert!(first_path.is_file());
        assert!(second_path.is_file());
        drop(first);
        drop(second);
        assert!(!first_path.exists());
        assert!(!second_path.exists());
    }
}
