#!/usr/bin/env python3
"""Golden-output comparison for `phyluce_assembly_match_contigs_to_probes`
(Python) vs. `phyluce assembly match-contigs-to-probes` (Rust).

Neither side actually invokes `lastz` here (it isn't installed in most CI/
sandbox environments): both implementations are pointed at the
precomputed `*.lastz` fixtures under
phyluce/tests/test-expected/probe-match/ via their `--skip-alignment`-
equivalent paths --

  * Python: the fixtures are copied into --output before running, and
    since `lastz.Align.run()` unconditionally re-runs lastz, this script
    instead calls the module's internals directly would be too invasive;
    so this harness only exercises the Rust `--skip-alignment` flag and
    verifies its result against the checked-in SQLite/CSV fixtures
    (which were themselves produced by a real run of the Python command).

This is a fixture-comparison harness (Rust output vs. known-good fixture),
not a live Python-vs-Rust subprocess diff like the other compare_*.py
scripts -- there is no lastz binary to run the legacy script against in
this environment.
"""
import shutil
import sqlite3
import subprocess
import sys
import tempfile
from pathlib import Path

from common import RUST_ROOT, find_fixture_repo
REPO_ROOT = find_fixture_repo()
FIXTURE_DIR = REPO_ROOT / "phyluce/tests/test-expected/probe-match"
CONTIGS_DIR = REPO_ROOT / "phyluce/tests/test-expected/spades/contigs"
PROBES = REPO_ROOT / "phyluce/tests/probes/uce-5k-probes.fasta"
EXPECTED_DB = FIXTURE_DIR / "probe.matches.sqlite"
EXPECTED_CSV = FIXTURE_DIR / "probe_match_results.csv"


def dump_db(path: Path):
    conn = sqlite3.connect(str(path))
    c = conn.cursor()
    cols = [
        row[1]
        for row in c.execute("PRAGMA table_info(matches)").fetchall()
        if row[1] != "uce"
    ]
    out = {}
    for table in ("matches", "match_map"):
        c.execute(f"SELECT uce,{','.join(cols)} FROM {table} ORDER BY uce")
        rows = c.fetchall()
        out[table] = {r[0]: dict(zip(cols, r[1:])) for r in rows}
    conn.close()
    return out


def read_csv_numeric_rows(path: Path):
    """organism -> [uce-contigs, total-contigs, dupe-probe-matches,
    loci-dropped, contigs-dropped], ignoring the header's first-column
    name ("taxon" in current source vs. "organism" in this older fixture)."""
    rows = {}
    with open(path) as f:
        next(f)  # header
        for line in f:
            parts = line.strip().split(",")
            rows[parts[0]] = parts[1:]
    return rows


def main():
    rust_bin = Path(
        sys.argv[1] if len(sys.argv) > 1 else RUST_ROOT / "target/debug/phyluce"
    )
    if not rust_bin.is_file():
        print(f"Rust binary not found at {rust_bin}", file=sys.stderr)
        return 2

    failed = 0
    with tempfile.TemporaryDirectory() as td:
        output_dir = Path(td) / "output"
        output_dir.mkdir()
        for lastz_file in FIXTURE_DIR.glob("*.lastz"):
            shutil.copy(lastz_file, output_dir / lastz_file.name)

        csv_out = Path(td) / "results.csv"
        proc = subprocess.run(
            [
                str(rust_bin), "assembly", "match-contigs-to-probes",
                "--contigs", str(CONTIGS_DIR),
                "--probes", str(PROBES),
                "--output", str(output_dir),
                "--skip-alignment",
                "--csv", str(csv_out),
            ],
            stdout=subprocess.PIPE, stderr=subprocess.PIPE, universal_newlines=True,
        )
        if proc.returncode != 0:
            print(f"Rust command failed: {proc.stdout}\n{proc.stderr}")
            return 1

        expected_db = dump_db(EXPECTED_DB)
        actual_db = dump_db(output_dir / "probe.matches.sqlite")
        for table in ("matches", "match_map"):
            keys_expected, keys_actual = set(expected_db[table]), set(actual_db[table])
            if keys_expected != keys_actual:
                failed += 1
                print(f"{table}: key set differs: {keys_expected ^ keys_actual}")
            for uce in keys_expected & keys_actual:
                if expected_db[table][uce] != actual_db[table][uce]:
                    failed += 1
                    print(
                        f"{table} {uce}: expected {expected_db[table][uce]} "
                        f"got {actual_db[table][uce]}"
                    )

        expected_csv = read_csv_numeric_rows(EXPECTED_CSV)
        actual_csv = read_csv_numeric_rows(csv_out)
        if expected_csv != actual_csv:
            failed += 1
            print(f"CSV rows differ: expected {expected_csv} got {actual_csv}")

    if failed:
        print(f"{failed} mismatch(es).")
        return 1
    print("match-contigs-to-probes: SQLite + CSV match the fixture exactly.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
