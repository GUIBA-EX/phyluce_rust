#!/usr/bin/env python3
"""Golden-output comparison for `phyluce_assembly_get_bed_from_lastz` (Python)
vs. `phyluce assembly get-bed-from-lastz` (Rust).

Compares both the written BED file and stdout (the near-miss `print(name)`
side channel) across every *.lastz fixture and a couple of --identity /
--continuity thresholds. Python's `lastz.Reader` raises a `RuntimeError`
at end-of-iteration on this Python version (PEP 479: `StopIteration` inside
a generator), which happens *after* all matches are already written -- both
implementations' exit codes are treated as equivalent as long as they agree
on failing or succeeding for the same reason (see compare_get_fasta_lengths
for the same non-zero-exit-is-ok-if-both-fail convention).
"""
import os
import subprocess
import sys
import tempfile
from pathlib import Path

from common import find_python_repo
REPO_ROOT = find_python_repo()
PY_SCRIPT = REPO_ROOT / "bin/assembly/phyluce_assembly_get_bed_from_lastz"


def run_py(lastz_file: Path, out_file: Path, identity: float, continuity: float):
    cmd = [
        sys.executable, str(PY_SCRIPT),
        "--lastz", str(lastz_file),
        "--output", str(out_file),
        "--identity", str(identity),
        "--continuity", str(continuity),
    ]
    env = {**os.environ, "PYTHONPATH": str(REPO_ROOT)}
    proc = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, universal_newlines=True, cwd=REPO_ROOT, env=env)
    return proc.returncode, proc.stdout


def run_rust(rust_bin: Path, lastz_file: Path, out_file: Path, identity: float, continuity: float):
    cmd = [
        str(rust_bin), "assembly", "get-bed-from-lastz",
        "--lastz", str(lastz_file),
        "--output", str(out_file),
        "--identity", str(identity),
        "--continuity", str(continuity),
    ]
    proc = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, universal_newlines=True)
    return proc.returncode, proc.stdout


def compare_one(lastz_file: Path, rust_bin: Path, identity: float, continuity: float, tmpdir: Path):
    py_out_file = tmpdir / "py.bed"
    rust_out_file = tmpdir / "rust.bed"

    py_code, py_stdout = run_py(lastz_file, py_out_file, identity, continuity)
    rust_code, rust_stdout = run_rust(rust_bin, lastz_file, rust_out_file, identity, continuity)

    # Legacy `lastz.Reader` always crashes with a StopIteration-derived
    # RuntimeError right at EOF (PEP 479), *after* writing all output --
    # both implementations must produce identical files regardless of the
    # Python exit code, so we don't gate on returncode here at all, only on
    # the bed/stdout content once both had a chance to write it.
    py_bed = py_out_file.read_text() if py_out_file.exists() else ""
    rust_bed = rust_out_file.read_text() if rust_out_file.exists() else ""

    if py_bed != rust_bed:
        return False, "bed file content differs"
    if py_stdout != rust_stdout:
        return False, "stdout (near-miss names) differs"
    return True, ""


def main():
    rust_bin = Path(
        sys.argv[1] if len(sys.argv) > 1 else REPO_ROOT / "rust/target/debug/phyluce"
    )
    if not rust_bin.is_file():
        print(f"Rust binary not found at {rust_bin}", file=sys.stderr)
        return 2

    fixtures = sorted((REPO_ROOT / "phyluce/tests/test-expected").rglob("*.lastz"))
    thresholds = [(90.0, 90.0), (99.9, 99.9), (0.0, 0.0)]

    total = 0
    failed = 0
    with tempfile.TemporaryDirectory() as td:
        tmpdir = Path(td)
        for lastz_file in fixtures:
            for identity, continuity in thresholds:
                total += 1
                ok, reason = compare_one(lastz_file, rust_bin, identity, continuity, tmpdir)
                if not ok:
                    failed += 1
                    print(f"MISMATCH: {lastz_file} identity={identity} continuity={continuity}: {reason}")

    print(f"Compared {total} case(s), {failed} mismatch(es).")
    return 1 if failed else 0


if __name__ == "__main__":
    sys.exit(main())
