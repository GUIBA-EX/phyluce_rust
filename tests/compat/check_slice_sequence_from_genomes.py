#!/usr/bin/env python3
"""Synthetic-regression check for `phyluce probe slice-sequence-from-genomes`.

Like `check_twobit_genome_sequences.py`, this can't be diffed against the
Python original: it needs `bx.seq.twobit`, which isn't installed here. This
hand-builds a tiny `.2bit` genome + a matching long-format lastz file and
checks the sliced-and-trimmed FASTA output against a hand-computed
expectation.
"""
import struct
import subprocess
import sys
import tempfile
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[3]
SIGNATURE = 0x1A412743
BASE_CODE = {"T": 0, "C": 1, "A": 2, "G": 3}


def build_twobit(name: str, bases: str) -> bytes:
    packed = bytearray()
    for i in range(0, len(bases), 4):
        chunk = bases[i : i + 4]
        byte = 0
        for j, ch in enumerate(chunk):
            byte |= BASE_CODE.get(ch.upper(), 0) << (6 - 2 * j)
        packed.append(byte)

    header = struct.pack("<IIII", SIGNATURE, 0, 1, 0)
    name_bytes = name.encode()
    index = struct.pack("<B", len(name_bytes)) + name_bytes
    seq_offset_pos = len(header) + len(index)
    index += struct.pack("<I", 0)
    record = struct.pack("<II", len(bases), 0)  # no N blocks
    record += struct.pack("<I", 0)  # no mask blocks
    record += struct.pack("<I", 0)  # reserved
    record += bytes(packed)

    seq_offset = len(header) + len(index)
    buf = bytearray(header + index)
    struct.pack_into("<I", buf, seq_offset_pos, seq_offset)
    return bytes(buf) + record


def lastz_long_line(name1, strand1, zstart1, end1, length1, name2, strand2, zstart2, end2, length2):
    diff = "." * 10
    return "\t".join(
        str(v)
        for v in [
            9000,
            f">{name1}",
            strand1,
            zstart1,
            end1,
            length1,
            f">{name2}",
            strand2,
            zstart2,
            end2,
            length2,
            diff,
            "10M",
            "10/10",
            "100.0%",
            "10/10",
            "100.0%",
            "10/10",
            "100.0%",
        ]
    ) + "\n"


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
        # 40bp contig: 5 Ns, 30bp ACGT-repeat, 5 Ns. The probe matches
        # [10, 20) within the clean region.
        bases = "N" * 5 + "ACGT" * 7 + "AA" + "N" * 5
        assert len(bases) == 40, len(bases)
        twobit_path = td / "genome.2bit"
        twobit_path.write_bytes(build_twobit("node1", bases))

        conf = td / "genomes.conf"
        conf.write_text(f"[chromos]\ntest1:{twobit_path}\n")

        lastz_dir = td / "lastz_in"
        lastz_dir.mkdir()
        (lastz_dir / "test1").write_text(
            lastz_long_line("node1", "+", 10, 20, 40, "uce-1_p1", "+", 0, 10, 10)
        )

        out_dir = td / "sliced_out"
        proc = subprocess.run(
            [
                str(rust_bin), "probe", "slice-sequence-from-genomes",
                "--conf", str(conf), "--lastz", str(lastz_dir), "--output", str(out_dir),
                "--flank", "5",
            ],
            capture_output=True, text=True,
        )
        out_fasta = out_dir / "test1.fasta"
        text = out_fasta.read_text() if proc.returncode == 0 and out_fasta.exists() else ""
        expected_slice = bases[5:25]  # ss=max(0,10-5)=5, se=min(40,20+5)=25
        if (
            proc.returncode != 0
            or "|contig:node1|slice:5-25|uce:uce-1|match:10-20|orient:+|probes:1" not in text
            or expected_slice not in text
        ):
            failed += 1
            print(
                f"slice-sequence-from-genomes: unexpected output {text!r} "
                f"stderr={proc.stderr!r}"
            )

    if failed:
        print(f"{failed} check(s) failed.")
        return 1
    print("slice-sequence-from-genomes (synthetic .2bit + lastz fixture): all checks passed.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
