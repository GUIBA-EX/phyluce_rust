#!/usr/bin/env python3
"""Golden-output comparison for `phyluce align
get-trimmed-alignments-from-untrimmed` (Rust) against the checked-in
`mafft-edge-trim/` fixture produced by
`phyluce_align_get_trimmed_alignments_from_untrimmed` (Python), driven from
the `mafft-for-edge-trim/` input fixture.

Pure internal algorithm (the native 3-stage phyluce trimming) -- no
external aligner binary required.
"""
import subprocess
import sys
import tempfile
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[3]
INPUT_DIR = REPO_ROOT / "phyluce/tests/test-expected/mafft-for-edge-trim"
EXPECTED_DIR = REPO_ROOT / "phyluce/tests/test-expected/mafft-edge-trim"


def main():
    rust_bin = Path(
        sys.argv[1] if len(sys.argv) > 1 else REPO_ROOT / "rust/target/debug/phyluce"
    )
    if not rust_bin.is_file():
        print(f"Rust binary not found at {rust_bin}", file=sys.stderr)
        return 2

    failed = 0
    with tempfile.TemporaryDirectory() as td:
        output_dir = Path(td) / "out"
        proc = subprocess.run(
            [
                str(rust_bin), "align", "get-trimmed-alignments-from-untrimmed",
                "--alignments", str(INPUT_DIR),
                "--output", str(output_dir),
            ],
            capture_output=True, text=True,
        )
        if proc.returncode != 0:
            print(f"command failed:\n{proc.stdout}\n{proc.stderr}")
            return 1

        expected_files = {p.name for p in EXPECTED_DIR.iterdir()}
        actual_files = {p.name for p in output_dir.iterdir()} if output_dir.exists() else set()
        if expected_files != actual_files:
            failed += 1
            print(f"file list differs: {expected_files ^ actual_files}")
        for name in expected_files & actual_files:
            if (EXPECTED_DIR / name).read_text() != (output_dir / name).read_text():
                failed += 1
                print(f"{name}: content differs")

    if failed:
        print(f"{failed} mismatch(es).")
        return 1
    print("get-trimmed-alignments-from-untrimmed: output matches mafft-edge-trim fixture exactly.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
