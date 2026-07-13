#!/usr/bin/env python3
"""Self-contained regression checks for the `phyluce genetrees` commands.

`dendropy` isn't installed in this environment, so there's no live Python
original to diff against here (unlike the other compare_*.py scripts).
Instead this locks in expected behavior against small synthetic trees:
rename-tree-leaves, get-tree-counts' rooting-invariant topology grouping,
get-mean-bootrep-support's mean/CI arithmetic, and the
generate+sort-multilocus-bootstrap-count round trip.
"""
import subprocess
import sys
import tempfile
from pathlib import Path

from common import find_python_repo
REPO_ROOT = find_python_repo()


def run(rust_bin, args, cwd=None):
    proc = subprocess.run([str(rust_bin), *args], capture_output=True, text=True, cwd=cwd)
    return proc.returncode, proc.stdout, proc.stderr


def main():
    rust_bin = Path(
        sys.argv[1] if len(sys.argv) > 1 else REPO_ROOT / "rust/target/debug/phyluce"
    ).resolve()
    if not rust_bin.is_file():
        print(f"Rust binary not found at {rust_bin}", file=sys.stderr)
        return 2

    failed = 0
    with tempfile.TemporaryDirectory() as td:
        td = Path(td)

        # --- rename-tree-leaves ---
        tree_in = td / "tree1.nwk"
        tree_in.write_text("(a,b,(c,d));\n")
        conf = td / "rename.conf"
        conf.write_text("[taxa]\na:Alligator_mississippiensis\nb:Gallus_gallus\n")
        tree_out = td / "renamed.nwk"
        code, _, err = run(rust_bin, [
            "genetrees", "rename-tree-leaves",
            "--input", str(tree_in), "--config", str(conf),
            "--output", str(tree_out), "--section", "taxa",
        ])
        expected = "(Alligator_mississippiensis,Gallus_gallus,(c,d));\n"
        if code != 0 or tree_out.read_text() != expected:
            failed += 1
            print(f"rename-tree-leaves: expected {expected!r}, got rc={code} err={err}")

        # --- get-tree-counts: same topology, different input rooting ---
        trees_dir = td / "treecounts"
        for name, newick in [
            ("locusA", "((a,b),c,d);"),
            ("locusB", "((c,d),a,b);"),
            ("locusC", "((a,c),b,d);"),
        ]:
            d = trees_dir / name
            d.mkdir(parents=True)
            (d / "RAxML_bestTree.FINAL").write_text(newick + "\n")
        support_out = td / "support.txt"
        code, stdout, err = run(rust_bin, [
            "genetrees", "get-tree-counts",
            "--trees", str(trees_dir), "--locus-support-output", str(support_out),
            "--root", "d", "--extension", "FINAL",
        ])
        if code != 0 or "2\t" not in stdout or "1\t" not in stdout:
            failed += 1
            print(f"get-tree-counts: unexpected output rc={code} stdout={stdout!r} err={err}")
        support_text = support_out.read_text()
        if "locusA" not in support_text or "locusB" not in support_text or "locusC" not in support_text:
            failed += 1
            print("get-tree-counts: locus-support-output missing expected loci")

        # --- get-mean-bootrep-support ---
        bootrep_dir = td / "bootrep"
        (bootrep_dir / "locusA").mkdir(parents=True)
        (bootrep_dir / "locusB").mkdir(parents=True)
        (bootrep_dir / "locusA" / "RAxML_bipartitions.FINAL").write_text(
            "((a:0.1,b:0.1)90:0.1,c:0.1,d:0.1);\n"
        )
        (bootrep_dir / "locusB" / "RAxML_bipartitions.FINAL").write_text(
            "((a:0.1,b:0.1)80:0.1,c:0.1,d:0.1);\n"
        )
        bootrep_conf = td / "bootrep.conf"
        bootrep_conf.write_text("[set1]\nlocusA\nlocusB\n")
        code, stdout, err = run(rust_bin, [
            "genetrees", "get-mean-bootrep-support",
            "--trees", str(bootrep_dir), "--config", str(bootrep_conf),
        ], cwd=td)
        if code != 0 or "set1,2,85,9.8" not in stdout:
            failed += 1
            print(f"get-mean-bootrep-support: expected mean 85/CI 9.8, got rc={code} stdout={stdout!r} err={err}")

        # --- generate + sort multilocus bootstrap round trip ---
        align_dir = td / "align"
        align_dir.mkdir()
        (align_dir / "uce-1.phylip").write_text("x")
        (align_dir / "uce-2.phylip").write_text("x")
        reps = td / "reps.txt"
        counts = td / "counts.txt"
        code, _, err = run(rust_bin, [
            "genetrees", "generate-multilocus-bootstrap-count",
            "--alignments", str(align_dir), "--bootstrap-replicates", str(reps),
            "--bootstrap-counts", str(counts), "--bootreps", "5",
        ])
        if code != 0 or not reps.exists() or not counts.exists():
            failed += 1
            print(f"generate-multilocus-bootstrap-count: command failed rc={code} err={err}")
        else:
            # counts.txt tells us exactly how many times random sampling
            # picked each locus -- generate at least that many synthetic
            # bootrep lines so the sort step never runs dry.
            needed = {}
            for line in counts.read_text().splitlines():
                path_str, count_str = line.rsplit(" ", 1)
                needed[Path(path_str).stem] = int(count_str)

            sortboot_in = td / "sortboot_in"
            for locus, count in needed.items():
                d = sortboot_in / locus
                d.mkdir(parents=True)
                (d / f"RAxML_bootstrap.{locus}.bootrep").write_text(
                    "\n".join(f"tree_{locus}_{i};" for i in range(count)) + "\n"
                )
            sortboot_out = td / "sortboot_out"
            code, _, err = run(rust_bin, [
                "genetrees", "sort-multilocus-bootstraps",
                "--input", str(sortboot_in), "--bootstrap-replicates", str(reps),
                "--output", str(sortboot_out),
            ])
            boot_files = sorted(sortboot_out.glob("boot*")) if sortboot_out.exists() else []
            if code != 0 or len(boot_files) != 5:
                failed += 1
                print(f"sort-multilocus-bootstraps: expected 5 boot files, rc={code} got={boot_files} err={err}")

    if failed:
        print(f"{failed} check(s) failed.")
        return 1
    print("genetrees (rename/tree-counts/bootrep-support/bootstrap round-trip): all synthetic checks passed.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
