#!/usr/bin/env python3
"""Golden-output comparison for the remaining P1 alignment commands (Rust)
against their fixtures:

- get-ry-recoded-alignments (plain + --binary)
- extract-taxa-from-alignments (--exclude + --include)
- split-concat-nexus-to-loci (round-trips our own concatenate-alignments
  output back to the original per-locus alignments)
- filter-alignments
- convert-one-align-to-another (fasta<->nexus, both directions)
"""
import subprocess
import sys
import tempfile
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[3]
E_DIR = REPO_ROOT / "phyluce/tests/test-expected"


def run(rust_bin, args):
    proc = subprocess.run([str(rust_bin), *args], capture_output=True, text=True)
    return proc.returncode, proc.stdout + proc.stderr


def compare_dir(expected_dir: Path, actual_dir: Path, label: str, failures: list):
    expected_files = {p.name for p in expected_dir.iterdir()}
    actual_files = {p.name for p in actual_dir.iterdir()} if actual_dir.exists() else set()
    if expected_files != actual_files:
        failures.append(f"{label}: file list differs: {expected_files ^ actual_files}")
        return
    for name in expected_files:
        if (expected_dir / name).read_text() != (actual_dir / name).read_text():
            failures.append(f"{label}: {name} content differs")


def main():
    rust_bin = Path(
        sys.argv[1] if len(sys.argv) > 1 else REPO_ROOT / "rust/target/debug/phyluce"
    )
    if not rust_bin.is_file():
        print(f"Rust binary not found at {rust_bin}", file=sys.stderr)
        return 2

    failures = []
    with tempfile.TemporaryDirectory() as td:
        td = Path(td)

        for label, extra, expected_name in [
            ("ry-recode", [], "mafft-gblocks-clean-ry"),
            ("ry-recode-binary", ["--binary"], "mafft-gblocks-clean-ry-binary"),
        ]:
            out = td / expected_name
            code, log = run(rust_bin, [
                "align", "get-ry-recoded-alignments",
                "--alignments", str(E_DIR / "mafft-gblocks-clean"),
                "--output", str(out),
                *extra,
            ])
            if code != 0:
                failures.append(f"{label}: command failed:\n{log}")
                continue
            compare_dir(E_DIR / expected_name, out, label, failures)

        exclude_out = td / "mafft-gblocks-clean-drop-gallus-gallus"
        code, log = run(rust_bin, [
            "align", "extract-taxa-from-alignments",
            "--alignments", str(E_DIR / "mafft-gblocks-clean"),
            "--output", str(exclude_out),
            "--input-format", "nexus", "--output-format", "nexus",
            "--exclude", "gallus_gallus",
        ])
        if code != 0:
            failures.append(f"extract-exclude: command failed:\n{log}")
        else:
            compare_dir(E_DIR / "mafft-gblocks-clean-drop-gallus-gallus", exclude_out, "extract-exclude", failures)

        include_out = td / "mafft-gblocks-clean-keep-gallus-and-peromyscus"
        code, log = run(rust_bin, [
            "align", "extract-taxa-from-alignments",
            "--alignments", str(E_DIR / "mafft-gblocks-clean"),
            "--output", str(include_out),
            "--input-format", "nexus", "--output-format", "nexus",
            "--include", "gallus_gallus", "peromyscus_maniculatus",
        ])
        if code != 0:
            failures.append(f"extract-include: command failed:\n{log}")
        else:
            compare_dir(E_DIR / "mafft-gblocks-clean-keep-gallus-and-peromyscus", include_out, "extract-include", failures)

        split_out = td / "split-back"
        code, log = run(rust_bin, [
            "align", "split-concat-nexus-to-loci",
            "--nexus", str(E_DIR / "mafft-gblocks-clean-concat-nexus/mafft-gblocks-clean-concat-nexus.nexus"),
            "--output", str(split_out),
        ])
        if code != 0:
            failures.append(f"split-concat: command failed:\n{log}")
        else:
            compare_dir(E_DIR / "mafft-gblocks-clean", split_out, "split-concat", failures)

        filter_out = td / "mafft-gblocks-filtered-alignments"
        code, log = run(rust_bin, [
            "align", "filter-alignments",
            "--alignments", str(E_DIR / "mafft-gblocks-clean"),
            "--output", str(filter_out),
            "--input-format", "nexus",
            "--containing-data-for", "gallus_gallus",
            "--min-length", "600", "--min-taxa", "3",
        ])
        if code != 0:
            failures.append(f"filter-alignments: command failed:\n{log}")
        else:
            expected_files = {
                p.name for p in (E_DIR / "mafft-gblocks-filtered-alignments").iterdir()
            }
            actual_files = {p.name for p in filter_out.iterdir()} if filter_out.exists() else set()
            if expected_files != actual_files:
                failures.append(f"filter-alignments: file list differs: {expected_files ^ actual_files}")

        f2n_out = td / "mafft-fasta-to-nexus"
        code, log = run(rust_bin, [
            "align", "convert-one-align-to-another",
            "--alignments", str(E_DIR / "mafft"),
            "--output", str(f2n_out),
            "--input-format", "fasta", "--output-format", "nexus",
        ])
        if code != 0:
            failures.append(f"convert-fasta-to-nexus: command failed:\n{log}")
        else:
            compare_dir(E_DIR / "mafft-fasta-to-nexus", f2n_out, "convert-fasta-to-nexus", failures)

        n2f_out = td / "mafft-nexus-to-fasta"
        code, log = run(rust_bin, [
            "align", "convert-one-align-to-another",
            "--alignments", str(f2n_out),
            "--output", str(n2f_out),
            "--input-format", "nexus", "--output-format", "fasta",
        ])
        if code != 0:
            failures.append(f"convert-nexus-to-fasta: command failed:\n{log}")
        else:
            compare_dir(E_DIR / "mafft-nexus-to-fasta", n2f_out, "convert-nexus-to-fasta", failures)

    if failures:
        for f in failures:
            print(f)
        print(f"{len(failures)} mismatch(es).")
        return 1
    print("ry-recode / extract-taxa / split-concat / filter-alignments / convert: all match fixtures exactly.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
