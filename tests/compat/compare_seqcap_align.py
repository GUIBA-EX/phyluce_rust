#!/usr/bin/env python3
"""Golden-output comparison for `phyluce align seqcap-align` (Rust, MAFFT
path) against the checked-in `mafft-no-trim/` fixture produced by
`phyluce_align_seqcap_align --aligner mafft --no-trim` (Python).

Requires an actual `mafft` binary (and `CONDA_PREFIX` pointing at its
install prefix, since `phyluce-config` resolves `[binaries] mafft` via
`$CONDA/bin/mafft`) -- both this script and the Rust CLI call the real
aligner. If mafft isn't available, this check is skipped (exit 0) rather
than failed, matching the rewrite plan's acknowledgment that external
binary availability varies by environment.
"""
import os
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path

from common import find_python_repo
REPO_ROOT = find_python_repo()
INPUT_FASTA = REPO_ROOT / "phyluce/tests/test-expected/taxon-set.incomplete.fasta"
EXPECTED_DIR = REPO_ROOT / "phyluce/tests/test-expected/mafft-no-trim"


def find_mafft_conda_prefix():
    mafft = shutil.which("mafft")
    if not mafft:
        return None
    # .../<conda-prefix>/bin/mafft -> <conda-prefix>
    return str(Path(mafft).resolve().parent.parent)


def main():
    rust_bin = Path(
        sys.argv[1] if len(sys.argv) > 1 else REPO_ROOT / "rust/target/debug/phyluce"
    )
    if not rust_bin.is_file():
        print(f"Rust binary not found at {rust_bin}", file=sys.stderr)
        return 2

    conda_prefix = os.environ.get("CONDA_PREFIX") or find_mafft_conda_prefix()
    if not conda_prefix or not shutil.which("mafft"):
        print("mafft not found on PATH; skipping seqcap-align golden comparison.")
        return 0

    env = {**os.environ, "CONDA_PREFIX": conda_prefix}
    failed = 0
    with tempfile.TemporaryDirectory() as td:
        output_dir = Path(td) / "out"
        proc = subprocess.run(
            [
                str(rust_bin), "align", "seqcap-align",
                "--input", str(INPUT_FASTA),
                "--output", str(output_dir),
                "--taxa", "4",
                "--no-trim",
            ],
            capture_output=True, text=True, env=env,
        )
        if proc.returncode != 0:
            print(f"command failed:\n{proc.stdout}\n{proc.stderr}")
            return 1

        expected_files = {p.name for p in EXPECTED_DIR.iterdir()}
        actual_files = {p.name for p in output_dir.iterdir()} if output_dir.exists() else set()
        if expected_files != actual_files:
            failed += 1
            print(f"file list differs: {expected_files ^ actual_files}")
        for name in expected_files & actual_files:
            if (EXPECTED_DIR / name).read_text() != (output_dir / name).read_text():
                failed += 1
                print(f"{name}: content differs")

    if failed:
        print(f"{failed} mismatch(es).")
        return 1
    print("seqcap-align (mafft, --no-trim): output matches mafft-no-trim fixture exactly.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
