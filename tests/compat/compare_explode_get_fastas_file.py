#!/usr/bin/env python3
"""Golden-output comparison for `phyluce assembly explode-get-fastas-file`
(Rust) against the checked-in `exploded-by-locus/` and `exploded-by-taxa/`
fixtures produced by `phyluce_assembly_explode_get_fastas_file` (Python).
"""
import subprocess
import sys
import tempfile
from pathlib import Path

from common import RUST_ROOT, find_fixture_repo
REPO_ROOT = find_fixture_repo()
INPUT_FASTA = REPO_ROOT / "phyluce/tests/test-expected/taxon-set.complete.fasta"
EXPECTED_BY_LOCUS = REPO_ROOT / "phyluce/tests/test-expected/exploded-by-locus"
EXPECTED_BY_TAXON = REPO_ROOT / "phyluce/tests/test-expected/exploded-by-taxa"


def run_rust(rust_bin, output_dir, by_taxon):
    cmd = [
        str(rust_bin), "assembly", "explode-get-fastas-file",
        "--input", str(INPUT_FASTA),
        "--output", str(output_dir),
    ]
    if by_taxon:
        cmd.append("--by-taxon")
    proc = subprocess.run(cmd, capture_output=True, text=True)
    return proc.returncode, proc.stdout + proc.stderr


def compare_dirs(expected_dir: Path, actual_dir: Path):
    mismatches = []
    expected_files = {p.name for p in expected_dir.iterdir()}
    actual_files = {p.name for p in actual_dir.iterdir()} if actual_dir.exists() else set()
    if expected_files != actual_files:
        mismatches.append(f"file list differs: {expected_files ^ actual_files}")
        return mismatches
    for name in expected_files:
        if (expected_dir / name).read_text() != (actual_dir / name).read_text():
            mismatches.append(f"{name}: content differs")
    return mismatches


def main():
    rust_bin = Path(
        sys.argv[1] if len(sys.argv) > 1 else RUST_ROOT / "target/debug/phyluce"
    )
    if not rust_bin.is_file():
        print(f"Rust binary not found at {rust_bin}", file=sys.stderr)
        return 2

    failed = 0
    with tempfile.TemporaryDirectory() as td:
        by_locus_out = Path(td) / "by-locus"
        code, log = run_rust(rust_bin, by_locus_out, by_taxon=False)
        if code != 0:
            failed += 1
            print(f"by-locus: command failed:\n{log}")
        else:
            for m in compare_dirs(EXPECTED_BY_LOCUS, by_locus_out):
                failed += 1
                print(f"by-locus: {m}")

        by_taxon_out = Path(td) / "by-taxon"
        code, log = run_rust(rust_bin, by_taxon_out, by_taxon=True)
        if code != 0:
            failed += 1
            print(f"by-taxon: command failed:\n{log}")
        else:
            for m in compare_dirs(EXPECTED_BY_TAXON, by_taxon_out):
                failed += 1
                print(f"by-taxon: {m}")

    if failed:
        print(f"{failed} mismatch(es).")
        return 1
    print("explode-get-fastas-file: by-locus + by-taxon outputs match fixtures exactly.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
