//! Validation for filenames derived from FASTA headers and configuration keys.

use std::path::{Component, Path, PathBuf};

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
}
