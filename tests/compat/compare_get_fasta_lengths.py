#!/usr/bin/env python3
"""Golden-output comparison: phyluce_assembly_get_fasta_lengths (Python)
vs. `phyluce assembly get-fasta-lengths` (Rust).

Per docs/rust-rewrite-plan.md section 4/11 ("Python/Rust golden output
comparison harness"). Runs both implementations over every *.fasta file
under phyluce/tests/test-expected/, in both human and --csv modes, and
reports mismatches.

Rules:
- If both implementations exit non-zero on malformed input, that's a pass
  regardless of the exact error text (both correctly reject bad input;
  Python's message is an implementation-detail traceback, Rust's is a
  clean diagnostic -- see docs/rust-rewrite-plan.md's compat notes).
- If exit codes differ, that's a failure.
- If both exit zero, integer fields must match exactly and float fields
  (avg/stderr/median) must match within a relative tolerance -- summation
  order differences (numpy pairwise summation vs. naive Rust summation)
  can differ in the last 1-2 ULPs without being a real bug.
"""
import math
import subprocess
import sys
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[3]
PY_SCRIPT = REPO_ROOT / "bin/assembly/phyluce_assembly_get_fasta_lengths"
FLOAT_TOL = 1e-6


def run_py(fasta: Path, csv: bool):
    cmd = [sys.executable, str(PY_SCRIPT), "--input", str(fasta)]
    if csv:
        cmd.append("--csv")
    proc = subprocess.run(
        cmd, capture_output=True, text=True, cwd=REPO_ROOT,
        env={**_base_env(), "PYTHONPATH": str(REPO_ROOT)},
    )
    return proc.returncode, proc.stdout + proc.stderr


def run_rust(rust_bin: Path, fasta: Path, csv: bool):
    cmd = [str(rust_bin), "assembly", "get-fasta-lengths", "--input", str(fasta)]
    if csv:
        cmd.append("--csv")
    proc = subprocess.run(cmd, capture_output=True, text=True)
    return proc.returncode, proc.stdout + proc.stderr


def _base_env():
    import os

    return dict(os.environ)


def parse_csv_row(row: str):
    return row.strip().split(",")


def fields_match(py_row, rust_row):
    if len(py_row) != len(rust_row):
        return False, "field count differs"
    for i, (p, r) in enumerate(zip(py_row, rust_row)):
        try:
            pf, rf = float(p), float(r)
            if not math.isclose(pf, rf, rel_tol=FLOAT_TOL, abs_tol=FLOAT_TOL):
                return False, f"field {i}: {p!r} vs {r!r}"
        except ValueError:
            if p != r:
                return False, f"field {i}: {p!r} vs {r!r}"
    return True, ""


def compare_one(fasta: Path, rust_bin: Path, csv: bool):
    py_code, py_out = run_py(fasta, csv)
    rust_code, rust_out = run_rust(rust_bin, fasta, csv)

    if py_code != 0 or rust_code != 0:
        if (py_code != 0) == (rust_code != 0):
            return True, ""
        return False, f"exit code mismatch: python={py_code} rust={rust_code}"

    if csv:
        ok, reason = fields_match(parse_csv_row(py_out), parse_csv_row(rust_out))
        return ok, reason

    # human report: compare line-by-line, tolerating float ULP differences
    py_lines = py_out.strip("\n").split("\n")
    rust_lines = rust_out.strip("\n").split("\n")
    if len(py_lines) != len(rust_lines):
        return False, "line count differs"
    for i, (pl, rl) in enumerate(zip(py_lines, rust_lines)):
        if pl == rl:
            continue
        p_val = pl.split("\t")[-1].replace(",", "")
        r_val = rl.split("\t")[-1].replace(",", "")
        ok, reason = fields_match([p_val], [r_val])
        if not ok:
            return False, f"line {i}: {pl!r} vs {rl!r} ({reason})"
    return True, ""


def main():
    rust_bin = Path(
        sys.argv[1] if len(sys.argv) > 1 else REPO_ROOT / "rust/target/debug/phyluce"
    )
    if not rust_bin.is_file():
        print(f"Rust binary not found at {rust_bin}; build with (cd rust && cargo build)", file=sys.stderr)
        return 2

    fixtures = sorted((REPO_ROOT / "phyluce/tests/test-expected").rglob("*.fasta"))
    total = 0
    failed = 0
    for fasta in fixtures:
        for csv in (False, True):
            total += 1
            ok, reason = compare_one(fasta, rust_bin, csv)
            if not ok:
                failed += 1
                print(f"MISMATCH: {fasta} csv={csv}: {reason}")

    print(f"Compared {total} case(s), {failed} mismatch(es).")
    return 1 if failed else 0


if __name__ == "__main__":
    sys.exit(main())
