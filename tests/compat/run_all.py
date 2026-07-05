#!/usr/bin/env python3
"""Run every Python/Rust golden-output comparison script in this directory.

Usage: rust/tests/compat/run_all.py [path-to-rust-binary]
"""
import subprocess
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent

SCRIPTS = [
    "compare_get_fasta_lengths.py",
    "compare_get_fastq_lengths.py",
    "compare_get_bed_from_lastz.py",
    "compare_get_bed_from_fasta.py",
    "compare_match_contigs_to_probes.py",
    "compare_get_match_counts.py",
    "compare_get_fastas_from_match_counts.py",
    "compare_explode_get_fastas_file.py",
    "compare_get_trimmed_alignments.py",
    "compare_seqcap_align.py",
    "compare_get_informative_sites.py",
    "compare_get_align_summary_data.py",
    "compare_concatenate_alignments.py",
    "compare_missing_data_and_remove_empty.py",
    "compare_p1_remaining.py",
    "compare_utilities_misc.py",
    "compare_chunk_fasta_for_ncbi.py",
    "check_genetrees_synthetic.py",
    "compare_probe_misc.py",
    "check_twobit_genome_sequences.py",
    "check_slice_sequence_from_genomes.py",
    "compare_align_assembly_batch2.py",
]


def main():
    rust_bin = sys.argv[1] if len(sys.argv) > 1 else None
    failed = []
    for script in SCRIPTS:
        cmd = [sys.executable, str(HERE / script)]
        if rust_bin:
            cmd.append(rust_bin)
        print(f"=== {script} ===")
        result = subprocess.run(cmd)
        if result.returncode != 0:
            failed.append(script)

    if failed:
        print(f"\nFAILED: {', '.join(failed)}")
        return 1
    print("\nAll golden-output comparisons passed.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
