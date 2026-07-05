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

The binary also supports the legacy script names for the second align/assembly
batch when it is invoked through a symlink or copied executable with the old
name. For example, a symlink named `phyluce_align_convert_degen_bases` is
translated internally to:

```text
phyluce align convert-degen-bases
```

Currently mapped legacy names:

```text
phyluce_align_convert_degen_bases
phyluce_align_explode_alignments
phyluce_align_extract_taxon_fasta_from_alignments
phyluce_align_format_concatenated_phylip_for_paml
phyluce_align_get_incomplete_matrix_estimates
phyluce_align_get_only_loci_with_min_taxa
phyluce_align_get_taxon_locus_counts_in_alignments
phyluce_align_move_align_by_conf_file
phyluce_align_randomly_sample_and_concatenate
phyluce_align_reduce_alignments_with_raxml
phyluce_align_remove_locus_name_from_files
phyluce_align_screen_alignments_for_problems
phyluce_align_get_smilogram_from_alignments
phyluce_assembly_screen_probes_for_dupes
phyluce_assembly_extract_contigs_to_barcodes
phyluce_assembly_match_contigs_to_barcodes
```

## Known differences

- `phyluce_assembly_match_contigs_to_barcodes`: the Rust port does not perform
  BOLD web lookups. Pass `--no-bold` to run the local LASTZ slicing step.
- Commands requiring external binaries such as LASTZ or RAxML still require the
  corresponding paths in `phyluce.conf`.
- Some Python 2 or Python 3-incompatible original scripts are matched to their
  intended behavior rather than their runtime failure mode.

## Logging

The Rust CLI accepts global logging flags on every command:

```text
--log-path DIR
--verbosity INFO|WARN|CRITICAL
```

When `--log-path` is supplied, the command writes `DIR/<program-name>.log`.
Native invocations write `phyluce.log`; legacy symlink invocations use the
legacy script name, for example `phyluce_align_convert_degen_bases.log`.

The log records command start/completion, version, raw arguments, and errors.
Normal command output remains on stdout/stderr so existing golden-output checks
and pipeline parsing are not affected. If `--log-path` is omitted, no log file
is created.

## Regression coverage

Compatibility checks live under `rust/tests/compat/`. The second align/assembly
batch uses checked-in `phyluce/tests/test-expected/` fixtures where available,
with synthetic smoke tests retained only for random, external-binary, or
Python-incompatible paths.
