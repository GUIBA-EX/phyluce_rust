//! Output-directory preparation and validation for data-derived filenames.

use std::path::{Component, Path, PathBuf};

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
}
