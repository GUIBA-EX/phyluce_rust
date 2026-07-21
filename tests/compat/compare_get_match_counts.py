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
import sqlite3
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
    proc = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, universal_newlines=True)
    return proc.returncode, proc.stdout + proc.stderr


def run_random_optimize(rust_bin: Path, database: Path, config: Path, output: Path):
    cmd = [
        str(rust_bin), "assembly", "get-match-counts",
        "--locus-db", str(database),
        "--taxon-list-config", str(config),
        "--taxon-group", "all",
        "--output", str(output),
        "--optimize", "--random",
        "--samples", "3",
        "--sample-size", "2",
        "--seed", "7",
    ]
    proc = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, universal_newlines=True)
    return proc.returncode, proc.stdout + proc.stderr


def write_optimize_fixture(root: Path):
    database = root / "optimize.sqlite"
    config = root / "optimize.conf"
    with sqlite3.connect(str(database)) as conn:
        conn.executescript(
            """
            CREATE TABLE matches (uce TEXT PRIMARY KEY, a INTEGER, b INTEGER, c INTEGER);
            CREATE TABLE match_map (uce TEXT PRIMARY KEY, a TEXT, b TEXT, c TEXT);
            INSERT INTO matches VALUES
                ('uce-1', 1, 1, 1),
                ('uce-2', 1, 1, 0),
                ('uce-3', 1, 1, 0),
                ('uce-4', 1, 0, 1);
            INSERT INTO match_map VALUES
                ('uce-1', 'a1(+)', 'b1(+)', 'c1(+)'),
                ('uce-2', 'a2(+)', 'b2(+)', ''),
                ('uce-3', 'a3(+)', 'b3(+)', ''),
                ('uce-4', 'a4(+)', '', 'c4(+)');
            """
        )
    config.write_text("[all]\na\nb\nc\n")
    return database, config


def main():
    rust_bin = Path(
        sys.argv[1] if len(sys.argv) > 1 else RUST_ROOT / "target/debug/phyluce"
    )
    if not rust_bin.is_file():
        print(f"Rust binary not found at {rust_bin}", file=sys.stderr)
        return 2

    failed = 0
    with tempfile.TemporaryDirectory() as td:
        temp_root = Path(td)
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

        database, config = write_optimize_fixture(temp_root)
        optimized = temp_root / "optimized.conf"
        code, log = run_random_optimize(rust_bin, database, config, optimized)
        expected_optimized = "[Organisms]\na\nb\n[Loci]\nuce-1\nuce-2\nuce-3\n"
        if code != 0:
            failed += 1
            print(f"random optimization failed:\n{log}")
        elif optimized.read_text() != expected_optimized:
            failed += 1
            print("random optimization output differs")
            print(f"  expected: {expected_optimized!r}")
            print(f"  actual:   {optimized.read_text()!r}")

    if failed:
        print(f"{failed} mismatch(es).")
        return 1
    print("get-match-counts: complete, incomplete, and random optimization outputs pass.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
