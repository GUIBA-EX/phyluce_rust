#!/usr/bin/env python3
"""Extraction parity for gzip input and detection of damaged trailing bytes."""
import gzip
import subprocess
import sys
import tempfile
from pathlib import Path

from common import LOCAL_FIXTURE_REPO, RUST_ROOT


def main():
    binary = Path(sys.argv[1] if len(sys.argv) > 1 else RUST_ROOT / "target/debug/phyluce").resolve()
    fixture = LOCAL_FIXTURE_REPO / "phyluce/tests/test-expected"
    with tempfile.TemporaryDirectory(prefix="phyluce-streaming-") as temp:
        root = Path(temp)
        contigs = root / "contigs"
        contigs.mkdir()
        sources = sorted((fixture / "spades/contigs").glob("*.fasta"))
        for source in sources:
            (contigs / source.name.replace(".contigs.fasta", ".fasta.gz")).write_bytes(gzip.compress(source.read_bytes()))

        def run(label, expect_success):
            output = root / (label + ".fasta")
            proc = subprocess.run([str(binary), "assembly", "get-fastas-from-match-counts",
                "--contigs", str(contigs), "--locus-db", str(fixture / "probe-match/probe.matches.sqlite"),
                "--match-count-output", str(fixture / "taxon-set.complete.conf"),
                "--output", str(output)], text=True, capture_output=True)
            assert (proc.returncode == 0) == expect_success, proc.stdout + proc.stderr
            if expect_success:
                assert output.read_bytes() == (fixture / "taxon-set.complete.fasta").read_bytes()

        run("gzip", True)
        target = contigs / sources[0].name.replace(".contigs.fasta", ".fasta.gz")
        # Invalid UTF-8 after every requested locus must still be checked.
        target.write_bytes(gzip.compress(sources[0].read_bytes() + b"\n>unused\n\xff\n"))
        run("invalid-utf8-tail", False)
        # A truncated trailer must not be ignored after selected records are emitted.
        target.write_bytes(gzip.compress(sources[0].read_bytes())[:-5])
        run("truncated-gzip", False)
    print("streaming extraction: gzip parity and trailing corruption checks pass.")


if __name__ == "__main__":
    main()
