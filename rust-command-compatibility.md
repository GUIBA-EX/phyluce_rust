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
- `prep-uce-align-files-for-ncbi` (`phyluce_ncbi_prep_uce_align_files_for_ncbi`)
  crashes on import against modern Biopython (`from Bio.Alphabet import IUPAC`
  -- `Bio.Alphabet` was removed). This port matches its intended behavior
  rather than that runtime failure.
- `get-match-counts --optimize` implements the intended exhaustive and random
  search behavior. Exhaustive rows are also written to the required output
  path, `--keep-counts` writes its CSV correctly, and Rust adds `--seed` for
  reproducible sampling.
- `sample-reads-from-files` shells out to `seqkit sample -p` (proportion)
  instead of the Python original's `seqtk sample -n` (exact count) --
  configure `[binaries] seqkit` in `phyluce.conf`. Same seed, different
  sampled reads (the two tools' algorithms/RNGs differ), and `-p` was a
  measured choice, not just a style pick: on a 604 MB/2M-read FASTQ,
  `seqkit -p` beat `seqtk -n` (3.1x faster, 1/5 the memory), while
  `seqkit -n` was both slower than `seqtk` and used ~15x its memory.
- `merge-multiple-gzip-files --trimmed` and `rename-tree-leaves --reroot` are
  implemented. `--reroot` matches DendroPy's `tree.reroot_at_node` semantics,
  including suppressing the unifurcations a naive reroot leaves behind.
  Several legacy alignment output formats are still not implemented and
  fail explicitly rather than silently changing behavior.
- `probe easy-stampy` replaces the hand-run `stampy.py` workflow from
  `docs/tutorials/tutorial-4.rst` with
  [probebwa](https://github.com/GUIBA-EX/probebwa), a stampy-compatible
  Rust mapper: `build-genome` → `build-hash` → `map` in one command, `--bam`
  for direct BAM output, and it skips a build step whose index file already
  exists (`--force-rebuild-index` to override). Configure the binary path
  under `[binaries] probebwa` in `phyluce.conf`. `probebwa` itself hasn't
  been validated on chromosome-scale (e.g. human) genomes yet -- see its own
  README for the E. coli-scale validation it has had.
- `match-contigs-to-probes`'s contig/probe name extraction has a hand-rolled
  fast path (`phyluce-assembly::fast_extract`, ~2.7x faster end to end) that
  only activates when `--regex`/`[headers]` are byte-identical to the
  packaged defaults; any customization keeps using the regex engine
  unchanged, and the fast path always falls back to it on a non-match rather
  than reporting one independently.
- `phyluce.conf`'s `[headers]` contig-naming patterns cover Trinity, Velvet,
  ABySS, IDBA, and SPAdes (matching the Python original) plus MEGAHIT
  (`k\d+_\d+`) and Flye (`contig_\d+`), which the original doesn't
  recognize. A contig header matching none of these patterns -- an
  assembler not on this list, or manually renamed contigs -- no longer
  aborts `match-contigs-to-probes`/`get-fastas-from-match-counts`: it falls
  back to the header's first whitespace-delimited token as the contig name
  and prints one aggregated warning per taxon (not per contig). Add a
  custom `[headers]` pattern if the fallback's guess is wrong for your data.

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
