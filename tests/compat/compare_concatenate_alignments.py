#!/usr/bin/env python3
"""Golden-output comparison for `phyluce align concatenate-alignments`
(Rust) against three fixture scenarios produced by
`phyluce_align_concatenate_alignments` (Python):
- NEXUS input -> NEXUS output (mafft-gblocks-clean -> mafft-gblocks-clean-concat-nexus)
- NEXUS input -> PHYLIP output (+ .charsets)
- FASTA input -> PHYLIP output (+ .charsets)
"""
import subprocess
import sys
import tempfile
from pathlib import Path

from common import find_python_repo
REPO_ROOT = find_python_repo()
E_DIR = REPO_ROOT / "phyluce/tests/test-expected"


def run_rust(rust_bin, alignments, input_format, output_dir, fmt_flag):
    cmd = [
        str(rust_bin), "align", "concatenate-alignments",
        "--alignments", str(alignments),
        "--input-format", input_format,
        "--output", str(output_dir),
        f"--{fmt_flag}",
    ]
    proc = subprocess.run(cmd, capture_output=True, text=True)
    return proc.returncode, proc.stdout + proc.stderr


def compare_dir(expected_dir: Path, actual_dir: Path):
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
        sys.argv[1] if len(sys.argv) > 1 else REPO_ROOT / "rust/target/debug/phyluce"
    )
    if not rust_bin.is_file():
        print(f"Rust binary not found at {rust_bin}", file=sys.stderr)
        return 2

    scenarios = [
        ("mafft-gblocks-clean", "nexus", "mafft-gblocks-clean-concat-nexus", "nexus"),
        ("mafft-gblocks-clean", "nexus", "mafft-gblocks-clean-concat-phylip", "phylip"),
        ("mafft-gblocks-clean-fasta", "fasta", "mafft-gblocks-clean-fasta-concat", "phylip"),
    ]

    failed = 0
    with tempfile.TemporaryDirectory() as td:
        for alignments_name, input_format, expected_name, fmt_flag in scenarios:
            alignments = E_DIR / alignments_name
            expected_dir = E_DIR / expected_name
            output_dir = Path(td) / expected_name
            code, log = run_rust(rust_bin, alignments, input_format, output_dir, fmt_flag)
            if code != 0:
                failed += 1
                print(f"{expected_name}: command failed:\n{log}")
                continue
            for m in compare_dir(expected_dir, output_dir):
                failed += 1
                print(f"{expected_name}: {m}")

    if failed:
        print(f"{failed} mismatch(es).")
        return 1
    print("concatenate-alignments: nexus/phylip/fasta-input scenarios all match fixtures exactly.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
