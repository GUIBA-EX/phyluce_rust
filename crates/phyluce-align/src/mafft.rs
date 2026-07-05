//! MAFFT wrapper mirroring `phyluce/mafft.py`'s `Align.run_alignment`.

use std::io::Write as _;

use phyluce_external::ExternalCommand;
use phyluce_io::read_fasta;

use crate::{Alignment, AlignmentRow};

#[derive(Debug, thiserror::Error)]
pub enum MafftError {
    #[error("{0}")]
    Io(#[from] std::io::Error),
    #[error("{0}")]
    External(#[from] phyluce_external::ExternalError),
    #[error("{0}")]
    Fasta(#[from] phyluce_io::FastaError),
}

/// Run MAFFT (`--adjustdirection --maxiterate 1000`) on the given
/// (id, sequence) records and return the resulting alignment. Mirrors
/// `Align.run_alignment`: MAFFT may prepend `_R_` to the id of any
/// sequence it reverse-complemented to fit the alignment -- this is
/// intentionally NOT stripped here (the legacy code doesn't either;
/// several downstream commands strip it explicitly when they need to).
pub fn run_mafft(mafft_bin: &str, records: &[(String, String)]) -> Result<Alignment, MafftError> {
    let dir = std::env::temp_dir();
    let input_path = dir.join(format!("phyluce-mafft-in-{}.fasta", std::process::id()));
    {
        let mut f = std::fs::File::create(&input_path)?;
        for (id, seq) in records {
            writeln!(f, ">{id}")?;
            writeln!(f, "{seq}")?;
        }
    }

    let report = ExternalCommand::new(mafft_bin)
        .args([
            "--adjustdirection".to_string(),
            "--maxiterate".to_string(),
            "1000".to_string(),
            input_path.to_string_lossy().to_string(),
        ])
        .run();
    let _ = std::fs::remove_file(&input_path);
    let report = report?;

    let out_path = dir.join(format!("phyluce-mafft-out-{}.fasta", std::process::id()));
    std::fs::write(&out_path, &report.stdout)?;
    let parsed = read_fasta(&out_path);
    let _ = std::fs::remove_file(&out_path);
    let parsed = parsed?;

    Ok(Alignment {
        rows: parsed
            .into_iter()
            .map(|r| AlignmentRow {
                id: r.id,
                seq: r.sequence.into_bytes(),
            })
            .collect(),
    })
}

pub fn mafft_binary_path(
    cfg: &phyluce_config::PhyluceConfig,
) -> Result<String, phyluce_config::ConfigError> {
    cfg.get_user_path("binaries", "mafft")
}
