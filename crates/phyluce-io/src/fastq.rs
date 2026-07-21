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
    #[error("{0}")]
    Utf8(#[from] std::str::Utf8Error),
}

/// Sequence-line lengths for one FASTQ (or FASTQ.gz) file: the length of the
/// first whitespace-delimited field on every 2nd line of each 4-line record
/// (header, sequence, plus, quality), matching `awk`'s `$1`.
///
/// Reads raw bytes via `read_until` into one reused buffer instead of
/// `BufRead::lines()`: `lines()` allocates a fresh `String` (and validates
/// it as UTF-8) for *every* line, including the 3-in-4 lines (header, `+`,
/// quality) this function never inspects. UTF-8 validation only happens for
/// the 1-in-4 sequence lines actually returned.
pub fn fastq_lengths(path: &Path) -> Result<Vec<usize>, FastqError> {
    let mut reader = open_maybe_gz(path)?;
    let mut lengths = Vec::new();
    let mut buf: Vec<u8> = Vec::new();
    let mut i = 0usize;
    loop {
        buf.clear();
        if reader.read_until(b'\n', &mut buf)? == 0 {
            break;
        }
        if i % 4 == 1 {
            let mut end = buf.len();
            while end > 0 && matches!(buf[end - 1], b'\n' | b'\r') {
                end -= 1;
            }
            let line = std::str::from_utf8(&buf[..end])?;
            let len = line.split_whitespace().next().unwrap_or("").len();
            lengths.push(len);
        }
        i += 1;
    }
    Ok(lengths)
}

/// Count FASTQ records without retaining their contents. This intentionally
/// mirrors the legacy line-count behaviour for incomplete trailing records.
///
/// Counts raw newlines via `read_until` into one reused `Vec<u8>` instead of
/// `BufRead::lines()`: since only the *count* of lines matters here, not
/// their content, this skips both the per-line `String` allocation and the
/// UTF-8 validation `lines()` would otherwise do for every line.
pub fn fastq_record_count(path: &Path) -> Result<usize, FastqError> {
    let mut reader = open_maybe_gz(path)?;
    let mut buf: Vec<u8> = Vec::new();
    let mut lines = 0usize;
    loop {
        buf.clear();
        if reader.read_until(b'\n', &mut buf)? == 0 {
            break;
        }
        lines += 1;
    }
    Ok(lines / 4)
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::fs::File;
    use std::io::Write;

    // Ad hoc benchmark for `fastq_lengths`'/`fastq_record_count`'s
    // read_until-based line reading. Run with:
    //   cargo +stable test --release -p phyluce-io --lib -- --ignored --nocapture bench_
    #[test]
    #[ignore]
    fn bench_fastq_lengths_large_file() {
        // ~200k reads * 150bp, resembling one lane of Illumina short-read
        // data -- the same shape used to compare against rust-bio (see the
        // note on `bench_read_fasta_large_file` in lib.rs).
        let dir = std::env::temp_dir().join("phyluce-io-bench");
        std::fs::create_dir_all(&dir).unwrap();
        let path = dir.join("bench.fastq");
        {
            let mut f = std::io::BufWriter::new(File::create(&path).unwrap());
            let bases = b"ACGT";
            for i in 0..200_000usize {
                writeln!(f, "@read_{i}/1").unwrap();
                let mut seq = [0u8; 150];
                for (j, b) in seq.iter_mut().enumerate() {
                    *b = bases[(i + j) % 4];
                }
                f.write_all(&seq).unwrap();
                f.write_all(b"\n+\n").unwrap();
                f.write_all(&[b'I'; 150]).unwrap();
                f.write_all(b"\n").unwrap();
            }
        }
        let file_len = std::fs::metadata(&path).unwrap().len();

        let start = std::time::Instant::now();
        let lengths = fastq_lengths(&path).unwrap();
        let elapsed = start.elapsed();
        eprintln!(
            "[bench] fastq_lengths: {} bytes / {} records in {:?} ({:.1} MB/s)",
            file_len,
            lengths.len(),
            elapsed,
            (file_len as f64 / 1_000_000.0) / elapsed.as_secs_f64()
        );

        let start = std::time::Instant::now();
        let count = fastq_record_count(&path).unwrap();
        let elapsed = start.elapsed();
        eprintln!(
            "[bench] fastq_record_count: {} bytes / {count} records in {:?} ({:.1} MB/s)",
            file_len,
            elapsed,
            (file_len as f64 / 1_000_000.0) / elapsed.as_secs_f64()
        );
        std::fs::remove_file(&path).ok();
    }

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
