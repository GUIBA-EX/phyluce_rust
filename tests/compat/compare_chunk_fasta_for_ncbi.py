#!/usr/bin/env python3
"""Golden-output comparison for `phyluce ncbi chunk-fasta-for-ncbi` (Rust)
against `phyluce_ncbi_chunk_fasta_for_ncbi` (Python), using a synthetic
25-record FASTA split into chunks of 10.
"""
import os
import subprocess
import sys
import tempfile
from pathlib import Path

from common import find_python_repo
REPO_ROOT = find_python_repo()


def main():
    rust_bin = Path(
        sys.argv[1] if len(sys.argv) > 1 else REPO_ROOT / "rust/target/debug/phyluce"
    ).resolve()
    if not rust_bin.is_file():
        print(f"Rust binary not found at {rust_bin}", file=sys.stderr)
        return 2

    with tempfile.TemporaryDirectory() as td:
        td = Path(td)
        input_fasta = td / "in.fasta"
        with open(input_fasta, "w") as f:
            for i in range(25):
                f.write(f">seq{i}\nACGT\n")

        py_dir = td / "py"
        py_dir.mkdir()
        env = {**os.environ, "PYTHONPATH": str(REPO_ROOT)}
        proc = subprocess.run(
            [
                sys.executable, str(REPO_ROOT / "bin/ncbi/phyluce_ncbi_chunk_fasta_for_ncbi"),
                "--input", str(input_fasta), "--chunk-size", "10", "--output-prefix", "split",
            ],
            capture_output=True, text=True, cwd=py_dir, env=env,
        )
        if proc.returncode != 0:
            print(f"Python command failed:\n{proc.stdout}\n{proc.stderr}")
            return 1

        rust_dir = td / "rust"
        rust_dir.mkdir()
        proc = subprocess.run(
            [
                str(rust_bin), "ncbi", "chunk-fasta-for-ncbi",
                "--input", str(input_fasta), "--chunk-size", "10", "--output-prefix", "split",
            ],
            capture_output=True, text=True, cwd=rust_dir,
        )
        if proc.returncode != 0:
            print(f"Rust command failed:\n{proc.stdout}\n{proc.stderr}")
            return 1

        py_files = sorted(p.name for p in py_dir.glob("split_*.fsa"))
        rust_files = sorted(p.name for p in rust_dir.glob("split_*.fsa"))
        if py_files != rust_files:
            print(f"file list differs: py={py_files} rust={rust_files}")
            return 1
        for name in py_files:
            if (py_dir / name).read_text() != (rust_dir / name).read_text():
                print(f"{name}: content differs")
                return 1

    print("chunk-fasta-for-ncbi: matches Python original exactly.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
