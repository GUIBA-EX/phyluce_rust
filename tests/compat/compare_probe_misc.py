#!/usr/bin/env python3
"""Golden-output / synthetic-regression checks for the newer
`phyluce probe` commands.

`remove-overlapping-probes-given-config` is diffed directly against the
Python original (no external tool dependency). The lastz-file-based
commands (`get-probe-bed-from-lastz-files`, `get-locus-bed-from-lastz-files`)
can't be diffed against Python here: `phyluce/lastz.py`'s `Reader` hits a
real PEP 479 bug at end-of-file (raises inside a generator), which for
these two specific commands means the Python original crashes *before*
writing any accumulated match data at all (see the bug report filed
alongside this port). So those two are checked against a synthetic
expected BED instead. The multi-fasta/multi-merge SQLite table commands
and get-subsets-of-tiled-probes are checked the same way (no fixtures
exist for them anywhere in the phyluce test suite either).
"""
import os
import shutil
import sqlite3
import subprocess
import sys
import tempfile
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[3]


def run_py(program, args):
    env = {**os.environ, "PYTHONPATH": str(REPO_ROOT)}
    cmd = [sys.executable, str(REPO_ROOT / "bin/probes" / program), *args]
    proc = subprocess.run(cmd, capture_output=True, text=True, env=env)
    return proc.returncode, proc.stdout


def run_rust(rust_bin, subcmd, args):
    proc = subprocess.run([str(rust_bin), "probe", subcmd, *args], capture_output=True, text=True)
    return proc.returncode, proc.stdout


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

        # --- remove-overlapping-probes-given-config: diffed against Python ---
        probes = td / "probes.fasta"
        probes.write_text(">uce-1_p1|source:x\nACGT\n>uce-2_p1|source:y\nACGT\n")
        conf = td / "exclude.conf"
        conf.write_text("[exclude]\nuce-2\n")
        py_out = td / "py_filtered.fasta"
        code, _ = run_py("phyluce_probe_remove_overlapping_probes_given_config", [
            "--probes", str(probes), "--config", str(conf), "--output", str(py_out),
        ])
        rust_out = td / "rust_filtered.fasta"
        rcode, _ = run_rust(rust_bin, "remove-overlapping-probes-given-config", [
            "--probes", str(probes), "--config", str(conf), "--output", str(rust_out),
        ])
        if code != 0 or rcode != 0 or py_out.read_text() != rust_out.read_text():
            failed += 1
            print("remove-overlapping-probes-given-config: mismatch")

        # --- get-probe-bed-from-lastz-files / get-locus-bed-from-lastz-files ---
        lastz_dir = td / "lastz_in"
        lastz_dir.mkdir()
        (lastz_dir / "something_v_alligator.lastz").write_text(
            "11415\t>NODE_6_length_856_cov_15.029963\t+\t363\t483\t120\t"
            ">uce-169_p2 |source:faircloth,probes-id:9855,probes-locus:169,probes-probe:2\t"
            "-\t0\t120\t120\t" + "." * 120 + "\t120M\t120/120\t100.0%\t120/120\t100.0%\t120/120\t100.0%\n"
            "11415\t>NODE_1_length_1016_cov_22.119667\t+\t516\t636\t120\t"
            ">uce-553_p1 |source:faircloth,probes-id:697,probes-locus:553,probes-probe:1\t"
            "-\t0\t120\t120\t" + "." * 120 + "\t120M\t120/120\t100.0%\t120/120\t100.0%\t120/120\t100.0%\n"
        )
        probe_bed_out = td / "probe_bed_out"
        rcode, _ = run_rust(rust_bin, "get-probe-bed-from-lastz-files", [
            "--alignments", str(lastz_dir), "--output", str(probe_bed_out),
        ])
        expected_probe_bed = (
            'track name="uce-v-alligator" description="UCE probe matches to alligator" visibility=2 itemRgb="On"\n'
            "NODE_6_length_856_cov_15.029963\t363\t483\tuce-169_p2\t1000\t+\t363\t483\t100,149,237\n"
            "NODE_1_length_1016_cov_22.119667\t516\t636\tuce-553_p1\t1000\t+\t516\t636\t100,149,237\n"
        )
        if rcode != 0 or (probe_bed_out / "alligator.probe.bed").read_text() != expected_probe_bed:
            failed += 1
            print("get-probe-bed-from-lastz-files: unexpected output")

        locus_bed_out = td / "locus_bed_out"
        rcode, _ = run_rust(rust_bin, "get-locus-bed-from-lastz-files", [
            "--alignments", str(lastz_dir), "--output", str(locus_bed_out),
        ])
        text = (locus_bed_out / "alligator.bed").read_text() if rcode == 0 else ""
        if rcode != 0 or "uce-169" not in text or "uce-553" not in text:
            failed += 1
            print("get-locus-bed-from-lastz-files: unexpected output")

        # --- get-multi-fasta-table + query-multi-fasta-table ---
        mft_dir = td / "mft_in"
        mft_dir.mkdir()
        (mft_dir / "taxonA.fasta").write_text(">uce-1|contig:x|coords:1-10|locus:1\nACGT\n")
        (mft_dir / "taxonB.fasta").write_text(
            ">uce-1|contig:y|coords:1-10|locus:1\nACGT\n>uce-2|contig:z|coords:1-10|locus:2\nACGT\n"
        )
        mft_db = td / "mft.sqlite"
        rcode, _ = run_rust(rust_bin, "get-multi-fasta-table", [
            "--fastas", str(mft_dir), "--output", str(mft_db), "--base-taxon", "taxonA",
        ])
        if rcode != 0:
            failed += 1
            print("get-multi-fasta-table: command failed")
        else:
            conn = sqlite3.connect(mft_db)
            rows = conn.execute("SELECT locus, taxonA, taxonB FROM taxonA ORDER BY locus").fetchall()
            conn.close()
            expected_rows = [("1", 1, 1), ("2", 0, 1)]
            if rows != expected_rows:
                failed += 1
                print(f"get-multi-fasta-table: unexpected rows {rows}")

        rcode, stdout = run_rust(rust_bin, "query-multi-fasta-table", [
            "--db", str(mft_db), "--base-taxon", "taxonA",
        ])
        if rcode != 0 or "Loci shared by 2 taxa:\t1" not in stdout:
            failed += 1
            print(f"query-multi-fasta-table: unexpected output {stdout!r}")

        # --- get-subsets-of-tiled-probes ---
        probes2 = td / "subset_probes.fasta"
        probes2.write_text(
            ">uce-1_p1 |source:x,probes-id:1,probes-locus:1,probes-probe:1,taxon:alligator\nACGT\n"
            ">uce-2_p1 |source:x,probes-id:2,probes-locus:2,probes-probe:1,taxon:gallus\nACGT\n"
        )
        subset_out = td / "subset_out.fasta"
        rcode, stdout = run_rust(rust_bin, "get-subsets-of-tiled-probes", [
            "--probes", str(probes2), "--taxa", "alligator", "--output", str(subset_out), "--regex", r"^(uce-\d+)(?:_p\d+.*)",
        ])
        if rcode != 0 or "uce-1_p1" not in subset_out.read_text() or "uce-2_p1" in subset_out.read_text():
            failed += 1
            print("get-subsets-of-tiled-probes: unexpected output")

        # --- get-screened-loci-by-proximity ---
        # locus 1 and 2 are within 10000bp of each other on chr1 (clustered,
        # only one survives); locus 3 is far away and always survives.
        proximity_in = td / "proximity_probes.fasta"
        proximity_in.write_text(
            ">uce-1_p1 |probes-global-chromo:chr1,probes-global-start:100,probes-global-end:200\nACGT\n"
            ">uce-2_p1 |probes-global-chromo:chr1,probes-global-start:300,probes-global-end:400\nACGT\n"
            ">uce-3_p1 |probes-global-chromo:chr1,probes-global-start:500000,probes-global-end:500100\nACGT\n"
        )
        proximity_out = td / "proximity_out.fasta"
        rcode, _ = run_rust(rust_bin, "get-screened-loci-by-proximity", [
            "--input", str(proximity_in), "--output", str(proximity_out), "--distance", "10000",
        ])
        out_text = proximity_out.read_text() if rcode == 0 else ""
        has_uce1 = ">uce-1_p1" in out_text
        has_uce2 = ">uce-2_p1" in out_text
        has_uce3 = ">uce-3_p1" in out_text
        if rcode != 0 or has_uce1 == has_uce2 or not has_uce3:
            failed += 1
            print(f"get-screened-loci-by-proximity: unexpected output {out_text!r}")

        # --- remove-duplicate-hits-from-probes-using-lastz ---
        # uce-1 self-to-self hits uce-2 as well (a cross-locus dupe hit) so
        # both should be dropped; uce-3 only hits itself so it survives.
        dupe_fasta = td / "dupe_probes.fasta"
        dupe_fasta.write_text(
            ">uce-1_p1\nACGT\n>uce-2_p1\nACGT\n>uce-3_p1\nACGT\n"
        )
        dupe_lastz = td / "dupe.lastz"
        dupe_lastz.write_text(
            "100\t>uce-1_p1\t+\t0\t4\t4\t>uce-1_p1\t+\t0\t4\t4\t"
            "....\t4M\t4/4\t100.0%\t4/4\t100.0%\t4/4\t100.0%\n"
            "100\t>uce-1_p1\t+\t0\t4\t4\t>uce-2_p1\t+\t0\t4\t4\t"
            "....\t4M\t4/4\t100.0%\t4/4\t100.0%\t4/4\t100.0%\n"
            "100\t>uce-3_p1\t+\t0\t4\t4\t>uce-3_p1\t+\t0\t4\t4\t"
            "....\t4M\t4/4\t100.0%\t4/4\t100.0%\t4/4\t100.0%\n"
        )
        rcode, _ = run_rust(rust_bin, "remove-duplicate-hits-from-probes-using-lastz", [
            "--fasta", str(dupe_fasta), "--lastz", str(dupe_lastz),
            "--probe-prefix", "uce-", "--long",
        ])
        dupe_out = td / "dupe_probes-DUPE-SCREENED.fasta"
        out_text = dupe_out.read_text() if rcode == 0 else ""
        if rcode != 0 or "uce-1_p1" in out_text or "uce-2_p1" in out_text or "uce-3_p1" not in out_text:
            failed += 1
            print(f"remove-duplicate-hits-from-probes-using-lastz: unexpected output {out_text!r}")

        # --- get-tiled-probe-from-multiple-inputs ---
        # A single 200bp locus from one taxon; sanity-check that probes are
        # designed at the expected length and global coordinates.
        tiled_dir = td / "tiled_in"
        tiled_dir.mkdir()
        locus_seq = "ACGT" * 50  # 200bp
        (tiled_dir / "taxonA.fasta").write_text(
            f">uce-1|contig:chr1|coords:1000-1200|locus:1\n{locus_seq}\n"
        )
        tiled_conf = td / "tiled_hits.conf"
        tiled_conf.write_text("[hits]\n1\n")
        tiled_out = td / "tiled_out.fasta"
        rcode, _ = run_rust(rust_bin, "get-tiled-probe-from-multiple-inputs", [
            "--fastas", str(tiled_dir), "--multi-fasta-output", str(tiled_conf),
            "--output", str(tiled_out), "--probe-prefix", "uce-",
            "--designer", "faircloth", "--design", "test-design",
            "--probe-length", "40", "--tiling-density", "2",
        ])
        out_text = tiled_out.read_text() if rcode == 0 else ""
        if rcode != 0 or "1_p1" not in out_text or "probes-global-chromo:chr1" not in out_text:
            failed += 1
            print(f"get-tiled-probe-from-multiple-inputs: unexpected output {out_text!r}")

        # --- get-tiled-probes ---
        tp_in = td / "tp_in.fasta"
        tp_in.write_text(f">uce-1|coords:1000-1200\n{locus_seq}\n")
        tp_out = td / "tp_out.fasta"
        tp_probe_bed = td / "tp_probes.bed"
        tp_locus_bed = td / "tp_loci.bed"
        rcode, _ = run_rust(rust_bin, "get-tiled-probes", [
            "--input", str(tp_in), "--output", str(tp_out),
            "--probe-prefix", "uce-", "--designer", "faircloth", "--design", "test-design",
            "--probe-length", "40", "--tiling-density", "2",
            "--probe-bed", str(tp_probe_bed), "--locus-bed", str(tp_locus_bed),
        ])
        out_text = tp_out.read_text() if rcode == 0 else ""
        pb_text = tp_probe_bed.read_text() if rcode == 0 else ""
        lb_text = tp_locus_bed.read_text() if rcode == 0 else ""
        if (
            rcode != 0
            or "uce-0_p1" not in out_text
            or not pb_text.startswith("track name=get_tiled_probes ")
            or not lb_text.startswith("track name=get_tiled_probes_loci ")
        ):
            failed += 1
            print(f"get-tiled-probes: unexpected output {out_text!r} {pb_text!r} {lb_text!r}")

        # --- reconstruct-uce-from-probe ---
        # Single-probe locus needs no aligner and should pass through
        # unchanged (minus the "_p1" suffix); a two-probe locus exercises
        # the MAFFT + dumb_consensus path if mafft is on PATH.
        recon_in = td / "recon_probes.fasta"
        recon_in.write_text(">uce-1_p1\nACGTACGTACGT\n")
        recon_out = td / "recon_out.fasta"
        rcode, _ = run_rust(rust_bin, "reconstruct-uce-from-probe", [
            "--input", str(recon_in), "--output", str(recon_out),
        ])
        out_text = recon_out.read_text() if rcode == 0 else ""
        if rcode != 0 or out_text.strip() != ">uce-1\nACGTACGTACGT":
            failed += 1
            print(f"reconstruct-uce-from-probe (single-probe): unexpected output {out_text!r}")

        if shutil.which("mafft"):
            recon_multi_in = td / "recon_multi_probes.fasta"
            recon_multi_in.write_text(">uce-2_p1\nACGTACGTACGTACGT\n>uce-2_p2\nACGTACGTACGTACGT\n")
            recon_multi_out = td / "recon_multi_out.fasta"
            rcode, _ = run_rust(rust_bin, "reconstruct-uce-from-probe", [
                "--input", str(recon_multi_in), "--output", str(recon_multi_out),
                "--mafft-binary", shutil.which("mafft"),
            ])
            multi_text = recon_multi_out.read_text() if rcode == 0 else ""
            if rcode != 0 or not multi_text.startswith(">uce-2\n") or "ACGT" not in multi_text:
                failed += 1
                print(f"reconstruct-uce-from-probe (multi-probe, mafft): unexpected output {multi_text!r}")

    if failed:
        print(f"{failed} check(s) failed.")
        return 1
    print("probe (remove-overlap/bed-from-lastz/multi-fasta-table/subsets): all checks passed.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
