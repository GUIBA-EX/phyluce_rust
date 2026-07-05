#!/usr/bin/env python3
"""Golden-output comparison for `phyluce_assembly_get_fastq_lengths` (Python)
vs. `phyluce assembly get-fastq-lengths` (Rust).

The legacy script's glob.glob() file order (and thus which file's basename
ends up in the CSV row) is filesystem-dependent/undefined; the Rust port
sorts files for determinism instead (see docs/rust-rewrite-plan.md section
9's list of quirks worth fixing). So this harness ignores the basename
field and only compares the numeric summary stats, with the same
ULP-tolerant float comparison as compare_get_fasta_lengths.py.
"""
import math
import os
import subprocess
import sys
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[3]
PY_SCRIPT = REPO_ROOT / "bin/assembly/phyluce_assembly_get_fastq_lengths"
FLOAT_TOL = 1e-6


def run_py(input_dir: Path, csv: bool):
    cmd = [sys.executable, str(PY_SCRIPT), "--input", str(input_dir)]
    if csv:
        cmd.append("--csv")
    env = {**os.environ, "PYTHONPATH": str(REPO_ROOT)}
    proc = subprocess.run(cmd, capture_output=True, text=True, cwd=REPO_ROOT, env=env)
    return proc.returncode, proc.stdout


def run_rust(rust_bin: Path, input_dir: Path, csv: bool):
    cmd = [str(rust_bin), "assembly", "get-fastq-lengths", "--input", str(input_dir)]
    if csv:
        cmd.append("--csv")
    proc = subprocess.run(cmd, capture_output=True, text=True)
    return proc.returncode, proc.stdout


def stats_fields_match(py_fields, rust_fields):
    # skip field 0 (basename) -- order is undefined upstream, fixed downstream
    if len(py_fields) != len(rust_fields):
        return False, "field count differs"
    for i in range(1, len(py_fields)):
        p, r = py_fields[i], rust_fields[i]
        try:
            if not math.isclose(float(p), float(r), rel_tol=FLOAT_TOL, abs_tol=FLOAT_TOL):
                return False, f"field {i}: {p!r} vs {r!r}"
        except ValueError:
            if p != r:
                return False, f"field {i}: {p!r} vs {r!r}"
    return True, ""


def compare_one(input_dir: Path, rust_bin: Path, csv: bool):
    py_code, py_out = run_py(input_dir, csv)
    rust_code, rust_out = run_rust(rust_bin, input_dir, csv)
    if py_code != 0 or rust_code != 0:
        if (py_code != 0) == (rust_code != 0):
            return True, ""
        return False, f"exit code mismatch: python={py_code} rust={rust_code}"

    if csv:
        # "All files in dir with <basename>,<count>,..." -- split off the
        # basename-bearing first field by comma same as the other fields.
        py_fields = py_out.strip().split(",")
        rust_fields = rust_out.strip().split(",")
        return stats_fields_match(py_fields, rust_fields)

    py_lines = py_out.strip("\n").split("\n")
    rust_lines = rust_out.strip("\n").split("\n")
    if len(py_lines) != len(rust_lines):
        return False, "line count differs"
    for i, (pl, rl) in enumerate(zip(py_lines, rust_lines)):
        if pl == rl:
            continue
        p_val = pl.split("\t")[-1].replace(",", "")
        r_val = rl.split("\t")[-1].replace(",", "")
        ok, reason = stats_fields_match(["_", p_val], ["_", r_val])
        if not ok:
            return False, f"line {i}: {pl!r} vs {rl!r} ({reason})"
    return True, ""


def main():
    rust_bin = Path(
        sys.argv[1] if len(sys.argv) > 1 else REPO_ROOT / "rust/target/debug/phyluce"
    )
    if not rust_bin.is_file():
        print(f"Rust binary not found at {rust_bin}", file=sys.stderr)
        return 2

    raw_reads_root = REPO_ROOT / "phyluce/tests/test-expected"
    fastq_dirs = sorted(
        {p.parent for p in raw_reads_root.rglob("*.fastq*")}
    )

    total = 0
    failed = 0
    for d in fastq_dirs:
        for csv in (False, True):
            total += 1
            ok, reason = compare_one(d, rust_bin, csv)
            if not ok:
                failed += 1
                print(f"MISMATCH: {d} csv={csv}: {reason}")

    print(f"Compared {total} case(s), {failed} mismatch(es).")
    return 1 if failed else 0


if __name__ == "__main__":
    sys.exit(main())
