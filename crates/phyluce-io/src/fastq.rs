//! FASTQ length extraction, mirroring the legacy
//! `phyluce_assembly_get_fastq_lengths`'s `gunzip | awk 'NR%4==2{print length($1)}'`
//! pipeline -- reimplemented directly instead of shelling out.

use std::io::{self, BufRead};
use std::path::Path;

use crate::open_maybe_gz;

#[derive(Debug, thiserror::Error)]
pub enum FastqError {
    #[error("{0}")]
    Io(#[from] io::Error),
}

/// Sequence-line lengths for one FASTQ (or FASTQ.gz) file: the length of the
/// first whitespace-delimited field on every 2nd line of each 4-line record
/// (header, sequence, plus, quality), matching `awk`'s `$1`.
pub fn fastq_lengths(path: &Path) -> Result<Vec<usize>, FastqError> {
    let reader = open_maybe_gz(path)?;
    let mut lengths = Vec::new();
    for (i, line) in reader.lines().enumerate() {
        let line = line?;
        if i % 4 == 1 {
            let len = line.split_whitespace().next().unwrap_or("").len();
            lengths.push(len);
        }
    }
    Ok(lengths)
}

/// Count FASTQ records without retaining their contents. This intentionally
/// mirrors the legacy line-count behaviour for incomplete trailing records.
pub fn fastq_record_count(path: &Path) -> Result<usize, FastqError> {
    let reader = open_maybe_gz(path)?;
    let mut lines = 0usize;
    for line in reader.lines() {
        line?;
        lines += 1;
    }
    Ok(lines / 4)
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::fs::File;
    use std::io::Write;

    #[test]
    fn extracts_sequence_line_lengths() {
        let dir = std::env::temp_dir().join("phyluce-io-fastq-tests");
        std::fs::create_dir_all(&dir).unwrap();
        let path = dir.join("reads.fastq");
        let mut f = File::create(&path).unwrap();
        f.write_all(b"@r1\nACGTACGT\n+\nIIIIIIII\n@r2\nACGT\n+\nIIII\n")
            .unwrap();
        assert_eq!(fastq_lengths(&path).unwrap(), vec![8, 4]);
        assert_eq!(fastq_record_count(&path).unwrap(), 2);
    }
}
