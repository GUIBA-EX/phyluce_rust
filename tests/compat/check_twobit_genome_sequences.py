#!/usr/bin/env python3
"""Synthetic-regression check for `phyluce probe get-genome-sequences-from-bed`.

`bx-python` (and thus `bx.seq.twobit`, which the Python original depends on)
isn't installed in this environment, so this can't be diffed against the
Python script directly. Instead we hand-encode a minimal, spec-correct
`.2bit` file in pure Python (no bx-python needed to *write* one) and check
the Rust CLI decodes/slices it correctly against a known-good expected
sequence.
"""
import struct
import subprocess
import sys
import tempfile
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[3]

SIGNATURE = 0x1A412743
BASE_CODE = {"T": 0, "C": 1, "A": 2, "G": 3}


def build_twobit(name: str, bases: str, n_blocks, mask_blocks) -> bytes:
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
    # placeholder offset patched below
    seq_offset_pos = len(header) + len(index)
    index += struct.pack("<I", 0)

    record = struct.pack("<II", len(bases), len(n_blocks))
    for s, _ in n_blocks:
        record += struct.pack("<I", s)
    for _, size in n_blocks:
        record += struct.pack("<I", size)
    record += struct.pack("<I", len(mask_blocks))
    for s, _ in mask_blocks:
        record += struct.pack("<I", s)
    for _, size in mask_blocks:
        record += struct.pack("<I", size)
    record += struct.pack("<I", 0)  # reserved
    record += bytes(packed)

    seq_offset = len(header) + len(index)
    buf = bytearray(header + index)
    struct.pack_into("<I", buf, seq_offset_pos, seq_offset)
    return bytes(buf) + record


def run_rust(rust_bin, subcmd, args):
    proc = subprocess.run([str(rust_bin), "probe", subcmd, *args], capture_output=True, text=True)
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
        # 20bp: 0-4 masked lowercase, 8-10 N-block, rest plain uppercase.
        bases = "acgtACGTNNACGTACGTAC"
        twobit_path = td / "genome.2bit"
        twobit_path.write_bytes(build_twobit("chr1", bases.upper(), [(8, 2)], [(0, 4)]))

        bed = td / "regions.bed"
        bed.write_text("chr1\t0\t20\n")
        out = td / "sliced.fasta"
        rcode, _, stderr = run_rust(rust_bin, "get-genome-sequences-from-bed", [
            "--bed", str(bed), "--twobit", str(twobit_path), "--output", str(out),
            "--filter-mask", "0.9", "--max-n", "5",
        ])
        text = out.read_text() if rcode == 0 else ""
        expected_seq = "acgtACGTNNACGTACGTAC"
        if rcode != 0 or expected_seq not in text or "chr1:0-20" not in text:
            failed += 1
            print(f"get-genome-sequences-from-bed: unexpected output {text!r} stderr={stderr!r}")

        # a stricter max-n of 0 should filter this region out entirely
        out2 = td / "sliced2.fasta"
        rcode, _, _ = run_rust(rust_bin, "get-genome-sequences-from-bed", [
            "--bed", str(bed), "--twobit", str(twobit_path), "--output", str(out2),
            "--max-n", "0",
        ])
        text2 = out2.read_text() if rcode == 0 else ""
        if rcode != 0 or text2.strip() != "":
            failed += 1
            print(f"get-genome-sequences-from-bed (--max-n 0 filter): unexpected output {text2!r}")

        # --- strip-masked-loci-from-set ---
        bed2 = td / "regions2.bed"
        bed2.write_text("chr1\t0\t20\nchr1\t0\t4\n")  # 2nd region is all masked
        strip_out = td / "stripped.bed"
        rcode, _, stderr = run_rust(rust_bin, "strip-masked-loci-from-set", [
            "--bed", str(bed2), "--twobit", str(twobit_path), "--output", str(strip_out),
            "--filter-mask", "0.5", "--max-n", "5",
        ])
        strip_text = strip_out.read_text() if rcode == 0 else ""
        if rcode != 0 or "chr1\t0\t20" not in strip_text or "chr1\t0\t4" in strip_text:
            failed += 1
            print(f"strip-masked-loci-from-set: unexpected output {strip_text!r} stderr={stderr!r}")

    if failed:
        print(f"{failed} check(s) failed.")
        return 1
    print("get-genome-sequences-from-bed (synthetic .2bit fixture): all checks passed.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
