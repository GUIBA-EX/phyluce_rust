#!/usr/bin/env python3
"""Self-contained CLI regression checks for alignment taxon renaming."""
import subprocess
import sys
import tempfile
from pathlib import Path

BINARY = Path(sys.argv[1] if len(sys.argv) > 1 else "target/debug/phyluce").resolve()


def main():
    with tempfile.TemporaryDirectory(prefix="phyluce-convert-names-") as temp:
        root = Path(temp)
        inputs = root / "input"
        inputs.mkdir()
        (inputs / "a.fasta").write_text(
            ">Homo_sapiens description\nACGT\n>Homo_sapient\nACGA\n"
        )
        (inputs / "b.fasta").write_text(
            ">Pan-troglodytes\nAGGT\n>Homo_sapiens\nACGT\n>abcdefghi\nAAAA\n"
        )

        def convert(output, *args, error=None):
            result = subprocess.run([
                str(BINARY), "align", "convert-one-align-to-another",
                "--alignments", str(inputs), "--output", str(root / output),
                "--input-format", "fasta", *args,
            ], text=True, capture_output=True)
            if error:
                assert result.returncode != 0, result
                assert error in result.stderr, result.stderr
                assert not (root / output).exists(), "invalid mappings created output"
            else:
                assert result.returncode == 0, result.stderr
            return root / output

        serial = convert("serial", "--shorten-names", "--cores", "1")
        parallel = convert("parallel", "--shorten-names", "--cores", "2")
        assert (serial / "a.fasta").read_text() == ">HomSap\nACGT\n>HomSap0\nACGA\n"
        assert (serial / "b.fasta").read_text() == ">PanTro\nAGGT\n>HomSap\nACGT\n>abcdef\nAAAA\n"
        for file in serial.iterdir():
            assert file.read_bytes() == (parallel / file.name).read_bytes()

        conf = root / "names.conf"
        valid = "[taxa]\nHomo_sapiens: HS\nHomo_sapient = HT\nPan-troglodytes: PT\nabcdefghi: ABC\n"
        conf.write_text(valid)
        mapped = convert("mapped", "--shorten-names", "--name-conf", str(conf), "--output-format", "nexus")
        # Round-trip through the actual NEXUS reader to check taxon identities and sequences.
        restored = root / "restored"
        subprocess.run([str(BINARY), "align", "convert-one-align-to-another",
                        "--alignments", str(mapped), "--output", str(restored)], check=True)
        assert (restored / "a.fasta").read_text() == ">HS\nACGT\n>HT\nACGA\n"
        assert (restored / "b.fasta").read_text() == ">PT\nAGGT\n>HS\nACGT\n>ABC\nAAAA\n"

        unchanged = convert("unchanged", "--name-conf", str(root / "nonexistent.conf"))
        assert ">Homo_sapiens\n" in (unchanged / "a.fasta").read_text()
        for index, (content, error) in enumerate([
            (valid.replace("Homo_sapiens", "homo_sapiens"), 'no name mapping for taxon "Homo_sapiens"'),
            (valid.replace("HT", "HS"), "duplicate short name"),
            (valid + "Homo_sapiens: XX\n", "duplicate taxon"),
            (valid.replace("[taxa]", "[names]"), "[taxa] section"),
            (valid.replace(": HS", ": "), "invalid short name"),
            (valid.replace(": HS", ": H S"), "invalid short name"),
        ]):
            conf.write_text(content)
            convert(f"invalid-{index}", "--shorten-names", "--name-conf", str(conf), error=error)
    print("convert-one-align-to-another: automatic and configured name mappings pass.")


if __name__ == "__main__":
    main()
