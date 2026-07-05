//! Multi-locus bootstrap replicate sampling, mirroring
//! `phyluce_genetrees_generate_multilocus_bootstrap_count` /
//! `phyluce_genetrees_sort_multilocus_bootstraps`.
//!
//! The legacy Python pair communicates via a `pickle`-serialized
//! `list[list[str]]`. Reproducing Python's pickle wire format byte-for-byte
//! isn't warranted here (it's a private hand-off between these two
//! phyluce commands, not a golden or user-facing file format) -- this
//! Rust port uses a plain, documented text format instead (one replicate
//! per line, loci comma-separated). **This means the Rust and Python
//! commands cannot be mixed** (generate with one, sort with the other);
//! use the same implementation for both ends of the pipeline.

use std::path::Path;

pub type Replicates = Vec<Vec<String>>;

pub fn write_replicates(path: &Path, replicates: &Replicates) -> std::io::Result<()> {
    let mut out = String::new();
    for replicate in replicates {
        out.push_str(&replicate.join(","));
        out.push('\n');
    }
    std::fs::write(path, out)
}

pub fn read_replicates(path: &Path) -> std::io::Result<Replicates> {
    let text = std::fs::read_to_string(path)?;
    Ok(text
        .lines()
        .filter(|l| !l.is_empty())
        .map(|l| l.split(',').map(str::to_string).collect())
        .collect())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn round_trips_replicates() {
        let dir = std::env::temp_dir().join("phyluce-genetrees-bootstrap-tests");
        std::fs::create_dir_all(&dir).unwrap();
        let path = dir.join(format!("reps-{}.txt", std::process::id()));
        let replicates = vec![
            vec!["a".to_string(), "b".to_string()],
            vec!["b".to_string(), "a".to_string()],
        ];
        write_replicates(&path, &replicates).unwrap();
        let read_back = read_replicates(&path).unwrap();
        assert_eq!(read_back, replicates);
    }
}
