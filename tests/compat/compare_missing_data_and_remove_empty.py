#!/usr/bin/env python3
"""Golden-output comparison for `phyluce align add-missing-data-designators`
and `phyluce align remove-empty-taxa` (Rust), chained the same way the
Python test suite does:

  mafft/*.fasta --add-missing-data-designators--> mafft-missing-data-designators/*.nexus
                --remove-empty-taxa-->             mafft-missing-data-designators-removed/*.nexus
"""
import subprocess
import sys
import tempfile
from pathlib import Path

from common import find_python_repo
REPO_ROOT = find_python_repo()
E_DIR = REPO_ROOT / "phyluce/tests/test-expected"


def compare_dir(expected_dir: Path, actual_dir: Path, label: str):
    mismatches = []
    expected_files = {p.name for p in expected_dir.iterdir()}
    actual_files = {p.name for p in actual_dir.iterdir()} if actual_dir.exists() else set()
    if expected_files != actual_files:
        mismatches.append(f"{label}: file list differs: {expected_files ^ actual_files}")
        return mismatches
    for name in expected_files:
        if (expected_dir / name).read_text() != (actual_dir / name).read_text():
            mismatches.append(f"{label}: {name} content differs")
    return mismatches


def main():
    rust_bin = Path(
        sys.argv[1] if len(sys.argv) > 1 else REPO_ROOT / "rust/target/debug/phyluce"
    )
    if not rust_bin.is_file():
        print(f"Rust binary not found at {rust_bin}", file=sys.stderr)
        return 2

    failed = 0
    with tempfile.TemporaryDirectory() as td:
        step1_out = Path(td) / "mafft-missing-data-designators"
        proc = subprocess.run(
            [
                str(rust_bin), "align", "add-missing-data-designators",
                "--alignments", str(E_DIR / "mafft"),
                "--output", str(step1_out),
                "--input-format", "fasta",
                "--match-count-output", str(E_DIR / "taxon-set.incomplete.conf"),
                "--incomplete-matrix", str(E_DIR / "taxon-set.incomplete"),
            ],
            capture_output=True, text=True,
        )
        if proc.returncode != 0:
            print(f"add-missing-data-designators failed:\n{proc.stdout}\n{proc.stderr}")
            return 1
        for m in compare_dir(E_DIR / "mafft-missing-data-designators", step1_out, "add-missing-data-designators"):
            failed += 1
            print(m)

        step2_out = Path(td) / "mafft-missing-data-designators-removed"
        proc = subprocess.run(
            [
                str(rust_bin), "align", "remove-empty-taxa",
                "--alignments", str(step1_out),
                "--output", str(step2_out),
                "--input-format", "nexus",
                "--output-format", "nexus",
            ],
            capture_output=True, text=True,
        )
        if proc.returncode != 0:
            print(f"remove-empty-taxa failed:\n{proc.stdout}\n{proc.stderr}")
            return 1
        for m in compare_dir(E_DIR / "mafft-missing-data-designators-removed", step2_out, "remove-empty-taxa"):
            failed += 1
            print(m)

    if failed:
        print(f"{failed} mismatch(es).")
        return 1
    print("add-missing-data-designators + remove-empty-taxa: outputs match fixtures exactly.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
