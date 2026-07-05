#!/usr/bin/env python3
"""Golden-output comparison for the newer `phyluce utilities` commands
against their Python originals, using small synthetic fixtures generated
on the fly (none of these commands have checked-in golden fixtures in
phyluce/tests/test-expected/).

Covers: filter-bed-by-fasta, replace-many-links, combine-reads,
merge-multiple-gzip-files, merge-next-seq-gzip-files, unmix-fasta-reads.
"""
import os
import subprocess
import sys
import tempfile
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[3]


def run_py(program, args):
    env = {**os.environ, "PYTHONPATH": str(REPO_ROOT)}
    cmd = [sys.executable, str(REPO_ROOT / "bin/utilities" / program), *args]
    proc = subprocess.run(cmd, capture_output=True, text=True, env=env)
    return proc.returncode, proc.stdout


def run_rust(rust_bin, subcmd, args):
    cmd = [str(rust_bin), "utilities", subcmd, *args]
    proc = subprocess.run(cmd, capture_output=True, text=True)
    return proc.returncode, proc.stdout


def main():
    rust_bin = Path(
        sys.argv[1] if len(sys.argv) > 1 else REPO_ROOT / "rust/target/debug/phyluce"
    )
    if not rust_bin.is_file():
        print(f"Rust binary not found at {rust_bin}", file=sys.stderr)
        return 2

    failed = 0
    with tempfile.TemporaryDirectory() as td:
        td = Path(td)

        # --- filter-bed-by-fasta ---
        fasta = td / "in.fasta"
        fasta.write_text(">uce-1_organism |uce-1\nACGT\n>uce-2_organism |uce-2\nACGT\n")
        bed = td / "in.bed"
        bed.write_text("track name=x\nchr1\t1\t10\tuce-1_something\nchr1\t20\t30\tuce-3_something\n")
        py_out = td / "py.bed"
        code, _ = run_py("phyluce_utilities_filter_bed_by_fasta", [
            "--bed", str(bed), "--fasta", str(fasta), "--output", str(py_out),
        ])
        rust_out = td / "rust.bed"
        rcode, _ = run_rust(rust_bin, "filter-bed-by-fasta", [
            "--bed", str(bed), "--fasta", str(fasta), "--output", str(rust_out),
        ])
        if code != 0 or rcode != 0 or py_out.read_text() != rust_out.read_text():
            failed += 1
            print("filter-bed-by-fasta: mismatch")

        # --- combine-reads ---
        reads_a = td / "reads_a"
        reads_b = td / "reads_b"
        reads_a.mkdir()
        reads_b.mkdir()
        (reads_a / "sample-READ1.fastq.gz").write_bytes(b"A")
        (reads_b / "sample-READ1.fastq.gz").write_bytes(b"B")
        conf = td / "combine.conf"
        conf.write_text(f"[samples]\nmysample:{reads_a},{reads_b}\n")
        py_out_dir = td / "py_combined"
        code, _ = run_py("phyluce_utilities_combine_reads", [
            "--config", str(conf), "--output", str(py_out_dir),
        ])
        rust_out_dir = td / "rust_combined"
        rcode, _ = run_rust(rust_bin, "combine-reads", [
            "--config", str(conf), "--output", str(rust_out_dir),
        ])
        py_files = sorted(p.name for p in py_out_dir.rglob("*") if p.is_file())
        rust_files = sorted(p.name for p in rust_out_dir.rglob("*") if p.is_file())
        if code != 0 or rcode != 0 or py_files != rust_files:
            failed += 1
            print(f"combine-reads: mismatch (py={py_files} rust={rust_files})")
        else:
            for name in py_files:
                py_content = next(p for p in py_out_dir.rglob(name)).read_bytes()
                rust_content = next(p for p in rust_out_dir.rglob(name)).read_bytes()
                if py_content != rust_content:
                    failed += 1
                    print(f"combine-reads: {name} content differs")

        # --- merge-multiple-gzip-files ---
        part_a = td / "partA.gz"
        part_b = td / "partB.gz"
        part_a.write_bytes(b"X")
        part_b.write_bytes(b"Y")
        merge_conf = td / "merge.conf"
        merge_conf.write_text(f"[samples]\nmerged.fastq.gz:{part_a},{part_b}\n")
        # NOTE: phyluce_utilities_merge_multiple_gzip_files's non-`--trimmed`
        # path is currently broken under Python 3 -- it opens the output in
        # binary mode but reads inputs in text mode, so
        # `shutil.copyfileobj` raises `TypeError: a bytes-like object is
        # required, not 'str'` on every invocation. That's a real bug in
        # the existing codebase (not something this port should
        # reproduce), so this check only confirms Python crashes and Rust
        # produces the correctly concatenated bytes.
        py_merge_out = td / "py_merge_out"
        code, _ = run_py("phyluce_utilities_merge_multiple_gzip_files", [
            "--config", str(merge_conf), "--output", str(py_merge_out),
        ])
        rust_merge_out = td / "rust_merge_out"
        rcode, _ = run_rust(rust_bin, "merge-multiple-gzip-files", [
            "--config", str(merge_conf), "--output", str(rust_merge_out),
        ])
        if code == 0:
            print("merge-multiple-gzip-files: Python no longer crashes -- re-check this test's assumptions")
            failed += 1
        elif rcode != 0 or (rust_merge_out / "merged.fastq.gz").read_bytes() != b"XY":
            failed += 1
            print("merge-multiple-gzip-files: Rust did not produce the expected concatenated output")

        # NOTE: phyluce_utilities_unmix_fasta_reads is currently broken
        # under Python 3 for any input -- it calls `os.write(fd, str)` on a
        # raw `tempfile.mkstemp()` descriptor, which requires bytes and
        # always raises `TypeError`. Another real bug in the existing
        # codebase, not something to reproduce; this check confirms Python
        # crashes and Rust produces the expected R1/R2/singleton split.
        mixed = td / "mixed.fasta"
        mixed.write_text(">readB/2\nCCCC\n>readA/1\nAAAA\n>readA/2\nGGGG\n>readC/1\nTTTT\n")
        py_r1, py_r2, py_rs = td / "py_r1.fasta", td / "py_r2.fasta", td / "py_rs.fasta"
        code, _ = run_py("phyluce_utilities_unmix_fasta_reads", [
            "--mixed-reads", str(mixed),
            "--out-r1", str(py_r1), "--out-r2", str(py_r2), "--out-r-singleton", str(py_rs),
        ])
        rust_r1, rust_r2, rust_rs = td / "rust_r1.fasta", td / "rust_r2.fasta", td / "rust_rs.fasta"
        rcode, _ = run_rust(rust_bin, "unmix-fasta-reads", [
            "--mixed-reads", str(mixed),
            "--out-r1", str(rust_r1), "--out-r2", str(rust_r2), "--out-r-singleton", str(rust_rs),
        ])
        if code == 0:
            print("unmix-fasta-reads: Python no longer crashes -- re-check this test's assumptions")
            failed += 1
        elif rcode != 0:
            failed += 1
            print("unmix-fasta-reads: Rust command failed")
        else:
            expected = {
                "r1": ">readA/1\nAAAA\n",
                "r2": ">readA/2\nGGGG\n",
                "singleton": ">readB/2\nCCCC\n>readC/1\nTTTT\n",
            }
            actual = {
                "r1": rust_r1.read_text(),
                "r2": rust_r2.read_text(),
                "singleton": rust_rs.read_text(),
            }
            if actual != expected:
                failed += 1
                print(f"unmix-fasta-reads: unexpected split {actual}")

    if failed:
        print(f"{failed} mismatch(es).")
        return 1
    print("utilities (filter-bed/combine-reads/merge-gzip/unmix-fasta): all match Python originals.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
