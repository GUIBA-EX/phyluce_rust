#!/usr/bin/env python3
"""Golden-output comparison for `phyluce_utilities_get_bed_from_fasta`
(Python) vs. `phyluce utilities get-bed-from-fasta` (Rust).

No repo fixture uses this command's expected
`id|contig:NAME|coords:BEGIN-END|locus:LOCUS` header shape, so this harness
generates a small synthetic FASTA covering the common cases (bare header,
extra whitespace after the id, and a --locus-prefix) instead.
"""
import os
import subprocess
import sys
import tempfile
from pathlib import Path

from common import find_python_repo
REPO_ROOT = find_python_repo()
PY_SCRIPT = REPO_ROOT / "bin/utilities/phyluce_utilities_get_bed_from_fasta"

FASTA_BODY = (
    ">uce-1|contig:NODE_1_length_500|coords:100-200|locus:uce-1\n"
    "ACGTACGTACGT\n"
    ">uce-2|contig:NODE_2_length_900|coords:0-450|locus:uce-2\n"
    "ACGTACGTACGT\n"
    ">uce-3 extra description|contig:NODE_3|coords:12-34|locus:uce-3\n"
    "ACGT\n"
)


def run_py(input_file: Path, out_file: Path, prefix: str):
    cmd = [
        sys.executable, str(PY_SCRIPT),
        "--input", str(input_file),
        "--output", str(out_file),
        "--locus-prefix", prefix,
    ]
    env = {**os.environ, "PYTHONPATH": str(REPO_ROOT)}
    proc = subprocess.run(cmd, capture_output=True, text=True, cwd=REPO_ROOT, env=env)
    return proc.returncode, proc.stdout + proc.stderr


def run_rust(rust_bin: Path, input_file: Path, out_file: Path, prefix: str):
    cmd = [
        str(rust_bin), "utilities", "get-bed-from-fasta",
        "--input", str(input_file),
        "--output", str(out_file),
        "--locus-prefix", prefix,
    ]
    proc = subprocess.run(cmd, capture_output=True, text=True)
    return proc.returncode, proc.stdout + proc.stderr


def main():
    rust_bin = Path(
        sys.argv[1] if len(sys.argv) > 1 else REPO_ROOT / "rust/target/debug/phyluce"
    )
    if not rust_bin.is_file():
        print(f"Rust binary not found at {rust_bin}", file=sys.stderr)
        return 2

    total = 0
    failed = 0
    with tempfile.TemporaryDirectory() as td:
        tmpdir = Path(td)
        input_file = tmpdir / "in.fasta"
        input_file.write_text(FASTA_BODY)

        for prefix in ("", "pre-"):
            total += 1
            py_out = tmpdir / "py.bed"
            rust_out = tmpdir / "rust.bed"
            py_code, _ = run_py(input_file, py_out, prefix)
            rust_code, _ = run_rust(rust_bin, input_file, rust_out, prefix)
            if py_code != 0 or rust_code != 0:
                if (py_code != 0) != (rust_code != 0):
                    failed += 1
                    print(f"MISMATCH prefix={prefix!r}: exit code mismatch python={py_code} rust={rust_code}")
                continue
            py_bed = py_out.read_text() if py_out.exists() else ""
            rust_bed = rust_out.read_text() if rust_out.exists() else ""
            if py_bed != rust_bed:
                failed += 1
                print(f"MISMATCH prefix={prefix!r}: bed content differs")
                print(f"  python: {py_bed!r}")
                print(f"  rust:   {rust_bed!r}")

    print(f"Compared {total} case(s), {failed} mismatch(es).")
    return 1 if failed else 0


if __name__ == "__main__":
    sys.exit(main())
