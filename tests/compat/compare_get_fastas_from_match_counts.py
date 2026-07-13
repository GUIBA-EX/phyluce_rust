#!/usr/bin/env python3
"""Golden-output comparison for `phyluce assembly get-fastas-from-match-counts`
(Rust) against the checked-in `taxon-set.complete.fasta` /
`taxon-set.incomplete.fasta` (+ `taxon-set.incomplete`) fixtures produced by
`phyluce_assembly_get_fastas_from_match_counts` (Python).
"""
import subprocess
import sys
import tempfile
from pathlib import Path

from common import RUST_ROOT, find_fixture_repo
REPO_ROOT = find_fixture_repo()
CONTIGS = REPO_ROOT / "phyluce/tests/test-expected/spades/contigs"
LOCUS_DB = REPO_ROOT / "phyluce/tests/test-expected/probe-match/probe.matches.sqlite"
EXPECTED_DIR = REPO_ROOT / "phyluce/tests/test-expected"


def run_rust(rust_bin, match_count_output, output, incomplete_out=None):
    cmd = [
        str(rust_bin), "assembly", "get-fastas-from-match-counts",
        "--contigs", str(CONTIGS),
        "--locus-db", str(LOCUS_DB),
        "--match-count-output", str(match_count_output),
        "--output", str(output),
    ]
    if incomplete_out:
        cmd.extend(["--incomplete-matrix", str(incomplete_out)])
    proc = subprocess.run(cmd, capture_output=True, text=True)
    return proc.returncode, proc.stdout + proc.stderr


def main():
    rust_bin = Path(
        sys.argv[1] if len(sys.argv) > 1 else RUST_ROOT / "target/debug/phyluce"
    )
    if not rust_bin.is_file():
        print(f"Rust binary not found at {rust_bin}", file=sys.stderr)
        return 2

    failed = 0
    with tempfile.TemporaryDirectory() as td:
        # complete matrix
        complete_out = Path(td) / "complete.fasta"
        code, log = run_rust(rust_bin, EXPECTED_DIR / "taxon-set.complete.conf", complete_out)
        if code != 0:
            failed += 1
            print(f"complete matrix: command failed:\n{log}")
        else:
            expected = (EXPECTED_DIR / "taxon-set.complete.fasta").read_text()
            actual = complete_out.read_text()
            if expected != actual:
                failed += 1
                print("complete matrix: FASTA output differs from fixture")

        # incomplete matrix
        incomplete_out = Path(td) / "incomplete.fasta"
        missing_out = Path(td) / "incomplete.missing"
        code, log = run_rust(
            rust_bin, EXPECTED_DIR / "taxon-set.incomplete.conf", incomplete_out, missing_out
        )
        if code != 0:
            failed += 1
            print(f"incomplete matrix: command failed:\n{log}")
        else:
            expected_fasta = (EXPECTED_DIR / "taxon-set.incomplete.fasta").read_text()
            actual_fasta = incomplete_out.read_text()
            if expected_fasta != actual_fasta:
                failed += 1
                print("incomplete matrix: FASTA output differs from fixture")
            expected_missing = (EXPECTED_DIR / "taxon-set.incomplete").read_text()
            actual_missing = missing_out.read_text()
            if expected_missing != actual_missing:
                failed += 1
                print("incomplete matrix: missing-locus report differs from fixture")

    if failed:
        print(f"{failed} mismatch(es).")
        return 1
    print("get-fastas-from-match-counts: complete + incomplete outputs match fixtures exactly.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
