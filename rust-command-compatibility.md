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
- `merge-multiple-gzip-files --trimmed` and `rename-tree-leaves --reroot` are
  implemented. `--trimmed` reuses `phyluce-assembly::raw_reads`'s R1/R2/
  singleton file discovery and writes `<output>/<name>/
  split-adapter-quality-trimmed/`; `--reroot` reroots at a leaf's parent by
  inverting the ancestor chain edge-by-edge, matching DendroPy's
  `tree.reroot_at_node` semantics (see `phyluce-genetrees::newick`).
  Several legacy alignment output formats are still not implemented and
  fail explicitly rather than silently changing behavior.
- `probe easy-stampy` replaces the hand-run `stampy.py` workflow from
  `docs/tutorials/tutorial-4.rst` with
  [probebwa](https://github.com/GUIBA-EX/probebwa), a stampy-compatible
  Rust mapper, chaining `build-genome` → `build-hash` → `map`. `--bam`
  writes BAM directly (no manual `samtools view` step). Configure the
  binary path under `[binaries] probebwa` in `phyluce.conf`. `probebwa`
  itself hasn't been validated on chromosome-scale (e.g. human) genomes yet
  -- see its own README for the E. coli-scale validation it has had.
- `match-contigs-to-probes`'s contig/probe name extraction has a
  hand-rolled fast path (`phyluce-assembly::fast_extract`) that activates
  only when `--regex`/`[headers]` are byte-identical to the packaged
  defaults and the input is ASCII; any customization keeps using the
  general regex engine unchanged. Benchmarked at ~2.7x end to end (300k
  synthetic LASTZ rows: ~175ms -> ~63ms). The fast path only ever
  *confirms* a match -- a `None` always falls through to `Regex::captures`
  rather than independently reporting "no match" -- so a bug there can only
  cost performance, never produce a wrong name; verified with differential
  fuzz tests (`fast_extract_*_matches_regex_oracle`) across 8 seeds and
  valid/boundary/garbage/non-ASCII inputs. The same benchmark suite found
  that swapping `ahash` in for the standard library's `HashMap`/`HashSet`
  is not worth it on its own: hashing is only ~10% of the pipeline's time,
  so the end-to-end gain was within run-to-run noise -- regex matching, not
  hashing, was the actual bottleneck.

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
