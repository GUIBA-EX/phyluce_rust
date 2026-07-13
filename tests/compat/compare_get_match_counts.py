#!/usr/bin/env python3
"""Golden-output comparison for `phyluce assembly get-match-counts` (Rust)
against the checked-in fixtures produced by
`phyluce_assembly_get_match_counts` (Python), for both complete and
incomplete matrix modes.

Like compare_match_contigs_to_probes.py, this compares the Rust output
directly to known-good fixtures rather than shelling out to the legacy
script, since the fixtures were generated against a specific
`probe.matches.sqlite` snapshot rather than anything this harness
regenerates on the fly.
"""
import subprocess
import sys
import tempfile
from pathlib import Path

from common import RUST_ROOT, find_fixture_repo
REPO_ROOT = find_fixture_repo()
LOCUS_DB = REPO_ROOT / "phyluce/tests/test-expected/probe-match/probe.matches.sqlite"
TAXON_CONFIG = REPO_ROOT / "phyluce/tests/test-conf/taxon-set.conf"
EXPECTED_COMPLETE = REPO_ROOT / "phyluce/tests/test-expected/taxon-set.complete.conf"
EXPECTED_INCOMPLETE = REPO_ROOT / "phyluce/tests/test-expected/taxon-set.incomplete.conf"


def run_rust(rust_bin: Path, output: Path, incomplete: bool):
    cmd = [
        str(rust_bin), "assembly", "get-match-counts",
        "--locus-db", str(LOCUS_DB),
        "--taxon-list-config", str(TAXON_CONFIG),
        "--taxon-group", "all",
        "--output", str(output),
    ]
    if incomplete:
        cmd.append("--incomplete-matrix")
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
        for incomplete, expected in (
            (False, EXPECTED_COMPLETE),
            (True, EXPECTED_INCOMPLETE),
        ):
            out = Path(td) / f"out-{incomplete}.conf"
            code, log = run_rust(rust_bin, out, incomplete)
            if code != 0:
                failed += 1
                print(f"incomplete={incomplete}: Rust command failed:\n{log}")
                continue
            actual = out.read_text()
            expected_text = expected.read_text()
            if actual != expected_text:
                failed += 1
                print(f"incomplete={incomplete}: output differs from {expected}")
                print(f"  expected: {expected_text!r}")
                print(f"  actual:   {actual!r}")

    if failed:
        print(f"{failed} mismatch(es).")
        return 1
    print("get-match-counts: complete + incomplete matrix outputs match fixtures exactly.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
