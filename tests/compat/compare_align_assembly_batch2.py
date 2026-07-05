#!/usr/bin/env python3
"""Golden-output / synthetic-regression checks for the second batch of
`phyluce align`/`phyluce assembly` commands (convert-degen-bases,
explode-alignments, extract-taxon-fasta-from-alignments,
format-concatenated-phylip-for-paml, get-incomplete-matrix-estimates,
get-only-loci-with-min-taxa, get-taxon-locus-counts-in-alignments,
move-align-by-conf-file, remove-locus-name-from-files,
screen-alignments-for-problems, screen-probes-for-dupes,
extract-contigs-to-barcodes).

`randomly-sample-and-concatenate` and `get-smilogram-from-alignments`
aren't diffed against Python here (the former is randomness-driven, the
latter would need a much larger fixture to exercise meaningfully) --
both get a basic Rust-only smoke check instead.
`reduce-alignments-with-raxml` and `match-contigs-to-barcodes` need
`raxml`/`lastz`, neither installed here, so they're skipped entirely
(matching the rest of this port's handling of unavailable binaries).
"""
import os
import subprocess
import sys
import tempfile
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[3]
EXPECTED = REPO_ROOT / "phyluce/tests/test-expected"


def run_py(program, args, subdir="align"):
    env = {**os.environ, "PYTHONPATH": str(REPO_ROOT)}
    cmd = [sys.executable, str(REPO_ROOT / "bin" / subdir / program), *args]
    proc = subprocess.run(cmd, capture_output=True, text=True, env=env)
    return proc.returncode, proc.stdout, proc.stderr


def run_rust(rust_bin, domain, subcmd, args):
    proc = subprocess.run([str(rust_bin), domain, subcmd, *args], capture_output=True, text=True)
    return proc.returncode, proc.stdout, proc.stderr


def read_fasta_records(path):
    records = {}
    name = None
    seq = []
    for raw in path.read_text().splitlines():
        if raw.startswith(">"):
            if name is not None:
                records[name] = "".join(seq)
            name = raw[1:].split()[0]
            seq = []
        elif raw.strip():
            seq.append(raw.strip())
    if name is not None:
        records[name] = "".join(seq)
    return records


def compare_dir_text(observed, expected):
    observed_files = sorted(p.name for p in observed.iterdir() if p.is_file())
    expected_files = sorted(p.name for p in expected.iterdir() if p.is_file())
    if observed_files != expected_files:
        return False, f"file set mismatch observed={observed_files} expected={expected_files}"
    for name in observed_files:
        if (observed / name).read_text() != (expected / name).read_text():
            return False, f"text mismatch in {name}"
    return True, ""


def compare_dir_fasta_records(observed, expected):
    observed_files = sorted(p.name for p in observed.iterdir() if p.is_file())
    expected_files = sorted(p.name for p in expected.iterdir() if p.is_file())
    if observed_files != expected_files:
        return False, f"file set mismatch observed={observed_files} expected={expected_files}"
    for name in observed_files:
        if read_fasta_records(observed / name) != read_fasta_records(expected / name):
            return False, f"FASTA record mismatch in {name}"
    return True, ""


NEXUS_A = """#NEXUS
begin data;
dimensions ntax=2 nchar=8;
format datatype=dna missing=? gap=-;
matrix
taxonA ACGTRYSW
taxonB ACGTACGT
;
end;
"""


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

        # --- convert-degen-bases ---
        py_out = td / "py_degen"
        code, _, err = run_py("phyluce_align_convert_degen_bases", [
            "--alignments", str(EXPECTED / "mafft-degen-bases"), "--output", str(py_out),
            "--input-format", "fasta", "--output-format", "nexus",
        ])
        rust_out = td / "rust_degen"
        rcode, _, rerr = run_rust(rust_bin, "align", "convert-degen-bases", [
            "--alignments", str(EXPECTED / "mafft-degen-bases"), "--output", str(rust_out),
            "--input-format", "fasta", "--output-format", "nexus",
        ])
        ok, why = compare_dir_text(rust_out, EXPECTED / "mafft-degen-bases-converted") if rcode == 0 else (False, "")
        if code != 0 or rcode != 0 or not ok:
            failed += 1
            print(f"convert-degen-bases: mismatch {why} (py rc={code} err={err!r}, rust rc={rcode} err={rerr!r})")

        # --- explode-alignments ---
        explode_out = td / "rust_explode_by_taxon"
        rcode, _, rerr = run_rust(rust_bin, "align", "explode-alignments", [
            "--alignments", str(EXPECTED / "mafft"), "--output", str(explode_out),
            "--input-format", "fasta", "--by-taxon",
        ])
        ok, why = compare_dir_fasta_records(explode_out, EXPECTED / "mafft-exploded-by-taxon") if rcode == 0 else (False, "")
        if rcode != 0 or not ok:
            failed += 1
            print(f"explode-alignments: mismatch {why} (rust rc={rcode} err={rerr!r})")

        # --- get-taxon-locus-counts-in-alignments ---
        py_csv = td / "py_counts.csv"
        code, _, err = run_py("phyluce_align_get_taxon_locus_counts_in_alignments", [
            "--alignments", str(EXPECTED / "mafft-gblocks-clean"), "--input-format", "nexus", "--output", str(py_csv),
        ])
        rust_csv = td / "rust_counts.csv"
        rcode, _, rerr = run_rust(rust_bin, "align", "get-taxon-locus-counts-in-alignments", [
            "--alignments", str(EXPECTED / "mafft-gblocks-clean"), "--input-format", "nexus", "--output", str(rust_csv),
        ])
        expected_counts = EXPECTED / "mafft-gblocks-clean-taxon-counts.csv"
        rust_rows = set(rust_csv.read_text().splitlines()) if rcode == 0 else set()
        expected_rows = set(expected_counts.read_text().splitlines())
        if code != 0 or rcode != 0 or rust_rows != expected_rows:
            failed += 1
            print(f"get-taxon-locus-counts-in-alignments: mismatch expected={expected_rows} rust={rust_rows} err={err!r}/{rerr!r}")

        # --- legacy argv[0] compatibility ---
        legacy_bin = td / "phyluce_align_get_taxon_locus_counts_in_alignments"
        os.symlink(rust_bin, legacy_bin)
        legacy_csv = td / "legacy_counts.csv"
        proc = subprocess.run([
            str(legacy_bin),
            "--alignments", str(EXPECTED / "mafft-gblocks-clean"),
            "--input-format", "nexus",
            "--output", str(legacy_csv),
        ], capture_output=True, text=True)
        legacy_rows = set(legacy_csv.read_text().splitlines()) if legacy_csv.exists() else set()
        if proc.returncode != 0 or legacy_rows != expected_rows:
            failed += 1
            print(f"legacy command name: mismatch rows={legacy_rows} err={proc.stderr!r}")

        # --- tracing/file logging ---
        logs_dir = td / "logs"
        logged_csv = td / "logged_counts.csv"
        rcode, _, rerr = run_rust(rust_bin, "align", "get-taxon-locus-counts-in-alignments", [
            "--alignments", str(EXPECTED / "mafft-gblocks-clean"),
            "--input-format", "nexus",
            "--output", str(logged_csv),
            "--log-path", str(logs_dir),
        ])
        log_file = logs_dir / "phyluce.log"
        log_text = log_file.read_text() if log_file.exists() else ""
        if rcode != 0 or "Starting phyluce" not in log_text or "Completed phyluce" not in log_text:
            failed += 1
            print(f"tracing log: missing expected log content rc={rcode} err={rerr!r} log={log_text!r}")

        # --- extract-taxon-fasta-from-alignments ---
        py_extract = td / "py_extract.fasta"
        code, _, err = run_py("phyluce_align_extract_taxon_fasta_from_alignments", [
            "--alignments", str(EXPECTED / "mafft-gblocks-clean"), "--taxon", "gallus_gallus", "--output", str(py_extract),
            "--input-format", "nexus",
        ])
        rust_extract = td / "rust_extract.fasta"
        rcode, _, rerr = run_rust(rust_bin, "align", "extract-taxon-fasta-from-alignments", [
            "--alignments", str(EXPECTED / "mafft-gblocks-clean"), "--taxon", "gallus_gallus", "--output", str(rust_extract),
            "--input-format", "nexus",
        ])
        expected_extract = EXPECTED / "mafft-gblocks-clean-gallus.fasta"
        if (
            code != 0 or rcode != 0
            or read_fasta_records(rust_extract) != read_fasta_records(expected_extract)
        ):
            failed += 1
            print(f"extract-taxon-fasta-from-alignments: mismatch (py rc={code} {err!r}, rust rc={rcode} {rerr!r})")
            print(f"  py={py_extract.read_text() if py_extract.exists() else None!r} rust={rust_extract.read_text() if rust_extract.exists() else None!r}")

        # --- screen-alignments-for-problems ---
        py_screen_out = td / "py_screen_out"
        code, pystdout, err = run_py("phyluce_align_screen_alignments_for_problems", [
            "--alignments", str(EXPECTED / "mafft-gblocks-clean-problems"), "--output", str(py_screen_out), "--input-format", "nexus",
        ])
        rust_screen_out = td / "rust_screen_out"
        rcode, _, rerr = run_rust(rust_bin, "align", "screen-alignments-for-problems", [
            "--alignments", str(EXPECTED / "mafft-gblocks-clean-problems"), "--output", str(rust_screen_out), "--input-format", "nexus",
        ])
        py_kept = sorted(p.name for p in py_screen_out.glob("*")) if code == 0 else None
        rust_kept = sorted(p.name for p in rust_screen_out.glob("*")) if rcode == 0 else None
        expected_kept = sorted(p.name for p in (EXPECTED / "mafft-gblocks-clean-problems-screened").glob("*"))
        if code != 0 or rcode != 0 or py_kept != rust_kept or rust_kept != expected_kept:
            failed += 1
            print(f"screen-alignments-for-problems: mismatch py_kept={py_kept} rust_kept={rust_kept} err={err!r}/{rerr!r}")

        # --- get-only-loci-with-min-taxa ---
        min_taxa_out_py = td / "py_min_taxa"
        code, _, err = run_py("phyluce_align_get_only_loci_with_min_taxa", [
            "--alignments", str(EXPECTED / "mafft-gblocks-clean"), "--taxa", "4", "--output", str(min_taxa_out_py),
            "--percent", "0.75", "--input-format", "nexus",
        ])
        min_taxa_out_rust = td / "rust_min_taxa"
        rcode, _, rerr = run_rust(rust_bin, "align", "get-only-loci-with-min-taxa", [
            "--alignments", str(EXPECTED / "mafft-gblocks-clean"), "--taxa", "4", "--output", str(min_taxa_out_rust),
            "--percent", "0.75", "--input-format", "nexus",
        ])
        py_kept = sorted(p.name for p in min_taxa_out_py.glob("*")) if code == 0 else None
        rust_kept = sorted(p.name for p in min_taxa_out_rust.glob("*")) if rcode == 0 else None
        expected_min_taxa = sorted(p.name for p in (EXPECTED / "mafft-gblocks-clean-75p").glob("*"))
        if code != 0 or rcode != 0 or py_kept != rust_kept or rust_kept != expected_min_taxa:
            failed += 1
            print(f"get-only-loci-with-min-taxa: mismatch py_kept={py_kept} rust_kept={rust_kept} err={err!r}/{rerr!r}")

        # --- remove-locus-name-from-files ---
        remove_out = td / "rust_remove_locus"
        rcode, _, rerr = run_rust(rust_bin, "align", "remove-locus-name-from-files", [
            "--alignments", str(EXPECTED / "mafft-gblocks"), "--output", str(remove_out),
            "--input-format", "nexus", "--output-format", "nexus",
        ])
        ok, why = compare_dir_text(remove_out, EXPECTED / "mafft-gblocks-clean") if rcode == 0 else (False, "")
        if rcode != 0 or not ok:
            failed += 1
            print(f"remove-locus-name-from-files: mismatch {why} (rust rc={rcode} err={rerr!r})")

        # --- move-align-by-conf-file ---
        move_in = td / "move_in"
        move_in.mkdir()
        (move_in / "uce-1.nex").write_text("dummy")
        (move_in / "uce-2.nex").write_text("dummy")
        move_conf = td / "move.conf"
        move_conf.write_text("[keep]\nuce-1.nex\n")
        py_move_out = td / "py_move_out"
        code, _, err = run_py("phyluce_align_move_align_by_conf_file", [
            "--conf", str(move_conf), "--alignments", str(move_in), "--output", str(py_move_out),
            "--extension", "nex",
        ])
        rust_move_out = td / "rust_move_out"
        rcode, _, rerr = run_rust(rust_bin, "align", "move-align-by-conf-file", [
            "--conf", str(move_conf), "--alignments", str(move_in), "--output", str(rust_move_out),
            "--extension", "nex",
        ])
        py_kept = sorted(p.name for p in py_move_out.glob("*")) if code == 0 else None
        rust_kept = sorted(p.name for p in rust_move_out.glob("*")) if rcode == 0 else None
        if code != 0 or rcode != 0 or py_kept != rust_kept:
            failed += 1
            print(f"move-align-by-conf-file: mismatch py_kept={py_kept} rust_kept={rust_kept} err={err!r}/{rerr!r}")

        # --- get-incomplete-matrix-estimates ---
        import sqlite3
        db_path = td / "matches.sqlite"
        conn = sqlite3.connect(db_path)
        conn.execute("CREATE TABLE matches (uce text, taxonA text, taxonB text)")
        conn.execute("INSERT INTO matches VALUES ('uce-1', '1', '1')")
        conn.execute("INSERT INTO matches VALUES ('uce-2', '1', '0')")
        conn.commit()
        conn.close()
        code, pyout, err = run_py("phyluce_align_get_incomplete_matrix_estimates", [
            "--db", str(db_path), "--min", "0", "--max", "1", "--step", "0.5",
        ])
        rcode, rustout, rerr = run_rust(rust_bin, "align", "get-incomplete-matrix-estimates", [
            "--db", str(db_path), "--min", "0", "--max", "1", "--step", "0.5",
        ])
        if code != 0 or rcode != 0:
            failed += 1
            print(f"get-incomplete-matrix-estimates: command failed py rc={code} {err!r}, rust rc={rcode} {rerr!r}")
        # (output format is a running total table; just check both ran cleanly)

        # --- format-concatenated-phylip-for-paml ---
        # NOTE: the Python original crashes under Python >=3.9 with
        # `TypeError: 'FakeSecHead' object is not iterable` --
        # `configparser`'s internals switched from calling `.readline()`
        # to iterating the file object directly, and `FakeSecHead` (a
        # local shim class used to fake an INI `[section]` header for a
        # headerless partition file) only implements `.readline()`. This
        # is a real, pre-existing Python-3 incompatibility bug in the
        # original script, not something introduced by this port -- so
        # we can't diff against it; just smoke-check the Rust CLI runs.
        phylip_in = td / "concat.phylip"
        phylip_in.write_text("2 8\ntaxonA   ACGTACGT\ntaxonB   ACGTACGT\n")
        partition_conf = td / "partitions.txt"
        partition_conf.write_text("DNA, p1 = 1-4\nDNA, p2 = 5-8\n")
        rust_paml = td / "rust_paml.phy"
        rcode, _, rerr = run_rust(rust_bin, "align", "format-concatenated-phylip-for-paml", [
            "--phylip-alignment", str(phylip_in), "--config", str(partition_conf), "--output", str(rust_paml),
        ])
        expected = "2 4\ntaxonA  ACGT\ntaxonB  ACGT\n\n2 4\ntaxonA  ACGT\ntaxonB  ACGT\n\n"
        if rcode != 0 or rust_paml.read_text() != expected:
            failed += 1
            print(f"format-concatenated-phylip-for-paml: unexpected output rc={rcode} {rerr!r}")
            print(f"  rust={rust_paml.read_text() if rust_paml.exists() else None!r}")

        # --- extract-contigs-to-barcodes ---
        contigs_dir = td / "contigs"
        contigs_dir.mkdir()
        (contigs_dir / "taxon-a.fasta").write_text(">NODE_1\nACGTACGT\n>NODE_2\nTTTT\n")
        barcode_conf = td / "barcodes.conf"
        barcode_conf.write_text("[assemblies]\ntaxon_a.fasta:NODE_1\n")
        py_barcodes = td / "py_barcodes.fasta"
        code, _, err = run_py("phyluce_assembly_extract_contigs_to_barcodes", [
            "--contigs", str(contigs_dir), "--config", str(barcode_conf), "--output", str(py_barcodes),
        ], subdir="assembly")
        rust_barcodes = td / "rust_barcodes.fasta"
        rcode, _, rerr = run_rust(rust_bin, "assembly", "extract-contigs-to-barcodes", [
            "--contigs", str(contigs_dir), "--config", str(barcode_conf), "--output", str(rust_barcodes),
        ])
        expected_barcode = ">taxon_a|NODE_1\nACGTACGT\n"
        if (
            code != 0 or rcode != 0
            or py_barcodes.read_text() != rust_barcodes.read_text()
            or rust_barcodes.read_text() != expected_barcode
        ):
            failed += 1
            print(f"extract-contigs-to-barcodes: mismatch py rc={code} {err!r}, rust rc={rcode} {rerr!r}")
            print(f"  py={py_barcodes.read_text() if py_barcodes.exists() else None!r}")
            print(f"  rust={rust_barcodes.read_text() if rust_barcodes.exists() else None!r}")

        # --- screen-probes-for-dupes (Rust-only: Python original is Python-2-only
        # syntax and can't run under Python 3 at all -- see module docs) ---
        dupe_lastz = td / "dupe_short.lastz"
        dupe_lastz.write_text(
            "100\t>uce-1|design:x\t+\t0\t4\t4\t>uce-1|design:x\t+\t0\t4\t4\t....\t4M\t4/4\t100.0%\t4/4\t100.0%\n"
            "100\t>uce-1|design:x\t+\t0\t4\t4\t>uce-2|design:x\t+\t0\t4\t4\t....\t4M\t4/4\t100.0%\t4/4\t100.0%\n"
        )
        rcode, rustout, rerr = run_rust(rust_bin, "assembly", "screen-probes-for-dupes", ["--lastz", str(dupe_lastz)])
        if rcode != 0 or "uce-1" not in rustout or "uce-2" not in rustout:
            failed += 1
            print(f"screen-probes-for-dupes: unexpected output {rustout!r} err={rerr!r}")

        # --- randomly-sample-and-concatenate (Rust-only smoke check) ---
        sample_dir = td / "sample_in"
        sample_dir.mkdir()
        for i in range(3):
            (sample_dir / f"uce-{i}.nex").write_text(NEXUS_A)
        sample_out = td / "sample_out"
        rcode, _, rerr = run_rust(rust_bin, "align", "randomly-sample-and-concatenate", [
            "--alignments", str(sample_dir), "--output", str(sample_out),
            "--sample-size", "2", "--replicates", "1",
        ])
        produced = list(sample_out.glob("*.nex")) if rcode == 0 else []
        if rcode != 0 or len(produced) != 1:
            failed += 1
            print(f"randomly-sample-and-concatenate: unexpected output rc={rcode} files={produced} err={rerr!r}")

        # --- get-smilogram-from-alignments (Rust-only smoke check) ---
        smilogram_out_file = td / "smilogram.csv"
        smilogram_out_missing = td / "smilogram_missing.csv"
        smilogram_db = td / "smilogram.sqlite"
        aln_dir = td / "aln_in"
        aln_dir.mkdir()
        (aln_dir / "uce-1.nexus").write_text(NEXUS_A)
        rcode, _, rerr = run_rust(rust_bin, "align", "get-smilogram-from-alignments", [
            "--alignments", str(aln_dir), "--output-file", str(smilogram_out_file),
            "--output-missing", str(smilogram_out_missing), "--output-database", str(smilogram_db),
            "--input-format", "nexus",
        ])
        if rcode != 0 or not smilogram_db.exists():
            failed += 1
            print(f"get-smilogram-from-alignments: command failed rc={rcode} err={rerr!r}")

    if failed:
        print(f"{failed} check(s) failed.")
        return 1
    print("align/assembly batch 2 (convert-degen-bases through extract-contigs-to-barcodes): all checks passed.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
