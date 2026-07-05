#!/usr/bin/env python3
"""Golden-output comparison for `phyluce align get-informative-sites`
(Rust) against `mafft-gblocks-clean-informative-sites.csv`, produced by
`phyluce_align_get_informative_sites` (Python) over the `mafft-gblocks-clean/`
fixture.

Row order is legitimately unordered upstream (`glob.glob`'s OS-dependent
order), so this compares rows as a locus -> fields dict, matching the
Python test's own comparison approach.
"""
import subprocess
import sys
import tempfile
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[3]
ALIGNMENTS = REPO_ROOT / "phyluce/tests/test-expected/mafft-gblocks-clean"
EXPECTED = REPO_ROOT / "phyluce/tests/test-expected/mafft-gblocks-clean-informative-sites.csv"


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
                str(rust_bin), "align", "get-informative-sites",
                "--alignments", str(ALIGNMENTS),
                "--input-format", "nexus",
                "--output", str(out),
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
    print("get-informative-sites: matches fixture exactly.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
