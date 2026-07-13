#!/usr/bin/env python3
"""Run self-contained golden fixture checks that do not execute Python phyluce."""

import subprocess
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
SCRIPTS = [
    "compare_match_contigs_to_probes.py",
    "compare_get_match_counts.py",
    "compare_get_fastas_from_match_counts.py",
    "compare_explode_get_fastas_file.py",
    "check_twobit_genome_sequences.py",
    "check_slice_sequence_from_genomes.py",
]


def main() -> int:
    rust_bin = sys.argv[1] if len(sys.argv) > 1 else None
    failed = []
    for script in SCRIPTS:
        cmd = [sys.executable, str(HERE / script)]
        if rust_bin:
            cmd.append(rust_bin)
        result = subprocess.run(cmd)
        if result.returncode != 0:
            failed.append(script)
    if failed:
        print(f"FAILED: {', '.join(failed)}")
        return 1
    print("All self-contained fixture checks passed.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
