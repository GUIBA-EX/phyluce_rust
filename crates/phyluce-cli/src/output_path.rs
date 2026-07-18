//! Output-directory preparation and validation for data-derived filenames.

use std::io::Write;
use std::path::{Component, Path, PathBuf};
use std::sync::atomic::{AtomicU64, Ordering};

static ATOMIC_WRITE_COUNTER: AtomicU64 = AtomicU64::new(0);

/// Create an output directory, or accept it only when it is already empty.
/// Non-empty directories are rejected so stale files cannot survive a rerun.
pub fn prepare_output_dir(directory: &Path) -> anyhow::Result<()> {
    if directory.exists() {
        anyhow::ensure!(
            directory.is_dir(),
            "{} is not a directory",
            directory.display()
        );
        anyhow::ensure!(
            std::fs::read_dir(directory)?.next().is_none(),
            "output directory {} is not empty; choose another directory or remove its contents",
            directory.display()
        );
    } else {
        std::fs::create_dir_all(directory)?;
    }
    Ok(())
}

/// Join one logical output filename beneath `directory` without allowing
/// absolute paths, parent traversal, or nested directories from input data.
pub fn output_file(directory: &Path, filename: &str) -> anyhow::Result<PathBuf> {
    let path = Path::new(filename);
    anyhow::ensure!(!filename.is_empty(), "output filename cannot be empty");
    let mut components = path.components();
    anyhow::ensure!(
        matches!(components.next(), Some(Component::Normal(_))) && components.next().is_none(),
        "unsafe output filename {filename:?}"
    );
    Ok(directory.join(path))
}

/// Remove a fixed-name sidecar from a previous external-tool invocation.
pub fn remove_stale_file(path: &Path) -> std::io::Result<()> {
    match std::fs::remove_file(path) {
        Ok(()) => Ok(()),
        Err(error) if error.kind() == std::io::ErrorKind::NotFound => Ok(()),
        Err(error) => Err(error),
    }
}

/// Reject an output path that resolves to an input file. This prevents a
/// successful-looking command from truncating its own database or config.
pub fn ensure_output_not_input(output: &Path, inputs: &[&Path]) -> anyhow::Result<()> {
    let output = std::fs::canonicalize(output).unwrap_or_else(|_| output.to_path_buf());
    for input in inputs {
        let input = std::fs::canonicalize(input).unwrap_or_else(|_| (*input).to_path_buf());
        anyhow::ensure!(
            output != input,
            "output path {} must not overwrite input {}",
            output.display(),
            input.display()
        );
    }
    Ok(())
}

/// Replace a text output only after the complete replacement has been written
/// beside it. The temporary file stays on the same filesystem as the target.
pub fn write_atomic(path: &Path, contents: impl AsRef<[u8]>) -> anyhow::Result<()> {
    let parent = path
        .parent()
        .filter(|parent| !parent.as_os_str().is_empty())
        .unwrap_or_else(|| Path::new("."));
    let name = path
        .file_name()
        .ok_or_else(|| anyhow::anyhow!("output path {} has no filename", path.display()))?
        .to_string_lossy();
    let sequence = ATOMIC_WRITE_COUNTER.fetch_add(1, Ordering::Relaxed);
    let temporary = parent.join(format!(
        ".{name}.phyluce-{}-{sequence}.tmp",
        std::process::id()
    ));

    let mut file = std::fs::OpenOptions::new()
        .write(true)
        .create_new(true)
        .open(&temporary)?;
    let write_result = (|| -> std::io::Result<()> {
        file.write_all(contents.as_ref())?;
        file.sync_all()
    })();
    drop(file);
    if let Err(error) = write_result {
        let _ = std::fs::remove_file(&temporary);
        return Err(error.into());
    }
    if let Err(error) = std::fs::rename(&temporary, path) {
        let _ = std::fs::remove_file(&temporary);
        return Err(error.into());
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn accepts_simple_filenames() {
        assert_eq!(
            output_file(Path::new("out"), "sample.fasta").unwrap(),
            PathBuf::from("out/sample.fasta")
        );
    }

    #[test]
    fn rejects_path_traversal_and_absolute_paths() {
        assert!(output_file(Path::new("out"), "../outside").is_err());
        assert!(output_file(Path::new("out"), "/tmp/outside").is_err());
        assert!(output_file(Path::new("out"), "nested/file").is_err());
    }

    #[test]
    fn rejects_nonempty_output_directories() {
        let directory =
            std::env::temp_dir().join(format!("phyluce-output-path-{}", std::process::id()));
        let _ = std::fs::remove_dir_all(&directory);
        prepare_output_dir(&directory).unwrap();
        prepare_output_dir(&directory).unwrap();
        std::fs::write(directory.join("stale.nex"), "x").unwrap();
        assert!(prepare_output_dir(&directory).is_err());
        std::fs::remove_dir_all(directory).unwrap();
    }

    #[test]
    fn removes_stale_sidecars_and_accepts_missing_files() {
        let path =
            std::env::temp_dir().join(format!("phyluce-stale-sidecar-{}", std::process::id()));
        std::fs::write(&path, "old").unwrap();
        remove_stale_file(&path).unwrap();
        assert!(!path.exists());
        remove_stale_file(&path).unwrap();
    }

    #[test]
    fn rejects_output_aliases_and_replaces_files_atomically() {
        let path = std::env::temp_dir().join(format!(
            "phyluce-atomic-write-{}-{}",
            std::process::id(),
            ATOMIC_WRITE_COUNTER.fetch_add(1, Ordering::Relaxed)
        ));
        std::fs::write(&path, "old").unwrap();
        assert!(ensure_output_not_input(&path, &[&path]).is_err());
        write_atomic(&path, "new").unwrap();
        assert_eq!(std::fs::read_to_string(&path).unwrap(), "new");
        std::fs::remove_file(path).unwrap();
    }
}
