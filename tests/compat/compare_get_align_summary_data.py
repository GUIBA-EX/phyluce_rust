#!/usr/bin/env python3
"""Golden-output comparison for `phyluce align get-align-summary-data`
(Rust) against `gblocks-clean-align-summary.csv`, produced by
`phyluce_align_get_align_summary_data --output-stats` (Python) over the
`mafft-gblocks-clean/` fixture. Row order is unordered upstream (same
rationale as compare_get_informative_sites.py).
"""
import subprocess
import sys
import tempfile
from pathlib import Path

from common import find_python_repo
REPO_ROOT = find_python_repo()
ALIGNMENTS = REPO_ROOT / "phyluce/tests/test-expected/mafft-gblocks-clean"
EXPECTED = REPO_ROOT / "phyluce/tests/test-expected/gblocks-clean-align-summary.csv"


def load(path: Path):
    rows = {}
    with open(path) as f:
        next(f)
        for line in f:
            parts = line.strip().split(",")
            rows[parts[0]] = parts[1:]
    return rows


def main():
    rust_bin = Path(
        sys.argv[1] if len(sys.argv) > 1 else REPO_ROOT / "rust/target/debug/phyluce"
    )
    if not rust_bin.is_file():
        print(f"Rust binary not found at {rust_bin}", file=sys.stderr)
        return 2

    with tempfile.TemporaryDirectory() as td:
        out = Path(td) / "out.csv"
        proc = subprocess.run(
            [
                str(rust_bin), "align", "get-align-summary-data",
                "--alignments", str(ALIGNMENTS),
                "--input-format", "nexus",
                "--output-stats", str(out),
            ],
            capture_output=True, text=True,
        )
        if proc.returncode != 0:
            print(f"command failed:\n{proc.stdout}\n{proc.stderr}")
            return 1
        actual = load(out)
        expected = load(EXPECTED)

    if actual != expected:
        print(f"key diff: {set(actual) ^ set(expected)}")
        for k in actual:
            if actual[k] != expected.get(k):
                print(f"{k}: expected {expected.get(k)} got {actual[k]}")
        return 1
    print("get-align-summary-data: matches fixture exactly.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
