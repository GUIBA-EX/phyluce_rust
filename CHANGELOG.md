# Changelog

## 0.3.0

### Added

- Automatic and configured taxon-name shortening for `convert-one-align-to-another` with `--shorten-names` and `--name-conf`.
- Parallel probe matching with `--cores`, fingerprint-validated `--resume`, and `--lastz-binary` selection.
- Configurable concatenation cache budget with `--memory-limit-mb` (default 256 MiB; 0 selects the disk path).
- Reproducible performance benchmarks and output-equivalence checks for extraction, matching, MAFFT scheduling, and concatenation.

### Improved

- Stream FASTA extraction and buffer output; reuse prepared SQLite statements during match ingestion.
- Prioritize larger MAFFT jobs while preserving alignment parameters, within-locus sequence order, and ordered reporting.
- Pack each locus into a contiguous sequence buffer with integer taxon identifiers.
- Parse concatenation inputs once. Spill over-budget data into an automatically cleaned snapshot and read sequence fragments directly when the matrix does not fit in memory, eliminating dense temporary disk matrices and their random writes.
- Preserve taxon order, missing characters, charset coordinates, and NEXUS/PHYLIP output bytes across storage paths.

### Performance and limits

On synthetic warm-cache Linux benchmarks, the latest disk-path change reduced concatenation time by roughly 29% for standard/short-locus inputs and reduced temporary disk space from 184 to 80 MB, 91 to 32 MB, and 110 to 10 MB across standard, short-locus, and sparse cases. These are stage measurements, not end-to-end guarantees.

Disk fragment indices reside in memory: short-locus disk-path RSS increased from approximately 23 to 53 MiB. The 48 MiB in-memory-matrix path showed a small timing regression and higher RSS. The cache budget is not a process RSS limit. See `tests/performance/README.md` for raw measurements, tradeoffs, and reproduction commands.

### Validation

Workspace tests, formatting, Clippy, self-contained compatibility fixtures, and upstream concatenation fixtures pass. The latest disk-output benchmarks include 288 output-equivalence runs covering FASTA/NEXUS input, NEXUS/PHYLIP output, and multiple memory budgets, with temporary-file cleanup checks.
