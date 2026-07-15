# Rust command compatibility

The Rust CLI is exposed as one binary:

```text
phyluce <domain> <command> [options]
```

For example:

```text
phyluce align convert-degen-bases --alignments in --output out
```

## Legacy command names

The binary supports all 74 legacy executable names when it is invoked through
a symlink or copied executable with the old
name. For example, a symlink named `phyluce_align_convert_degen_bases` is
translated internally to:

```text
phyluce align convert-degen-bases
```

Examples of mapped legacy names:

```text
phyluce_align_convert_degen_bases
phyluce_assembly_match_contigs_to_probes
phyluce_probe_easy_lastz
phyluce_genetrees_get_tree_counts
phyluce_utilities_combine_reads
phyluce_ncbi_chunk_fasta_for_ncbi
phyluce_workflow
```

## Known differences

- `phyluce_assembly_match_contigs_to_barcodes`: the Rust port does not perform
  BOLD web lookups. Pass `--no-bold` to run the local LASTZ slicing step.
- Commands requiring external binaries such as LASTZ or RAxML still require the
  corresponding paths in `phyluce.conf`.
- `reconstruct-uce-from-probe` uses MAFFT by default for multi-probe loci;
  `--muscle-binary` explicitly selects the legacy MUSCLE 3/Clustal path.
- Some Python 2 or Python 3-incompatible original scripts are matched to their
  intended behavior rather than their runtime failure mode.
- `get-match-counts --optimize` implements the intended exhaustive and random
  search behavior. Exhaustive rows are also written to the required output
  path, `--keep-counts` writes its CSV correctly, and Rust adds `--seed` for
  reproducible sampling.
- `merge-multiple-gzip-files --trimmed`, `rename-tree-leaves --reroot`, and
  several legacy alignment output formats are not implemented.
  These options fail explicitly rather than silently changing behavior.

## Logging

The Rust CLI accepts global logging flags on every command:

```text
--log-path DIR
--verbosity INFO|WARN|CRITICAL
```

When `--log-path` is supplied, the command writes `DIR/<program-name>.log`.
Native invocations write `phyluce.log`; legacy symlink invocations use the
legacy script name, for example `phyluce_align_convert_degen_bases.log`.

The log records command start/completion, version, raw arguments, errors, and
normal command diagnostics. Normal command output remains on stdout/stderr, so
existing golden-output checks and pipeline parsing are not affected. If
`--log-path` is omitted, no log file is created.

## Regression coverage

Compatibility checks live under `tests/compat/`. `run_fixtures.py` is fully
self-contained and uses fixtures checked into this repository. Live
Python/Rust comparisons use the original source tree selected with
`PHYLUCE_PYTHON_REPO` and may also require external bioinformatics tools.
