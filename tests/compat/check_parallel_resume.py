#!/usr/bin/env python3
"""Exercise bounded sample execution and validated resume with a controlled LASTZ stand-in.
The stand-in copies real LASTZ fixtures; it does not validate LASTZ itself.
"""
import os
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path

from common import LOCAL_FIXTURE_REPO, RUST_ROOT
from compare_match_contigs_to_probes import dump_db


def main():
    binary = Path(sys.argv[1] if len(sys.argv) > 1 else RUST_ROOT / "target/debug/phyluce").resolve()
    fixtures = LOCAL_FIXTURE_REPO / "phyluce/tests"
    with tempfile.TemporaryDirectory(prefix="phyluce-parallel-resume-") as td:
        root = Path(td)
        contigs = root / "contigs"
        shutil.copytree(fixtures / "test-expected/spades/contigs", contigs)
        probes = root / "probes.fasta"
        shutil.copyfile(fixtures / "probes/uce-5k-probes.fasta", probes)
        tool = root / "lastz"
        events = root / "events"
        tool.write_text("#!" + sys.executable + "\n" + '''import os, shutil, sys, time
from pathlib import Path
name = Path(sys.argv[1].split("[")[0]).stem
with open(os.environ["TEST_LASTZ_EVENTS"], "a") as log:
    log.write(name + "\\n")
output = Path(next(a.split("=", 1)[1] for a in sys.argv if a.startswith("--output=")))
if name == os.environ.get("TEST_LASTZ_FAIL"):
    output.write_text("partial failed output")
    sys.exit(3)
time.sleep(0.02)
shutil.copyfile(Path(os.environ["TEST_LASTZ_FIXTURES"]) / (name + ".lastz"), output)
''')
        tool.chmod(0o755)
        config = root / "phyluce.conf"
        config.write_text((RUST_ROOT / "config/phyluce.conf").read_text() + f"\n[binaries]\nlastz:{tool}\n")
        env = dict(os.environ, PHYLUCE_CONFIG=str(config), TEST_LASTZ_EVENTS=str(events),
                   TEST_LASTZ_FIXTURES=str(fixtures / "test-expected/probe-match"))

        def run(label, *args, fail=False, extra_env=None):
            output = root / label
            proc = subprocess.run([str(binary), "assembly", "match-contigs-to-probes",
                "--contigs", str(contigs), "--probes", str(probes), "--output", str(output),
                "--csv", str(root / (label + ".csv")), "--keep-duplicates", str(root / (label + ".dupes")),
                *([] if "--skip-alignment" in args else ["--lastz-binary", str(tool)]),
                *args], text=True, capture_output=True, env=env | (extra_env or {}))
            assert (proc.returncode != 0) == fail, proc.stdout + proc.stderr
            return output, proc

        serial, _ = run("serial")
        parallel, _ = run("parallel", "--cores", "3")
        assert dump_db(serial / "probe.matches.sqlite") == dump_db(parallel / "probe.matches.sqlite")
        assert (root / "serial.csv").read_bytes() == (root / "parallel.csv").read_bytes()
        assert (root / "serial.dupes").read_bytes() == (root / "parallel.dupes").read_bytes()
        for file in serial.glob("*.lastz"):
            assert file.read_bytes() == (parallel / file.name).read_bytes()
        # Verify the real fixture output under parallel parsing too.
        run("parallel", "--skip-alignment", "--cores", "3")
        assert dump_db(serial / "probe.matches.sqlite") == dump_db(parallel / "probe.matches.sqlite")

        def count():
            return len(events.read_text().splitlines())

        cached, _ = run("cached", "--resume", "--cores", "2")
        before = count()
        _, proc = run("cached", "--resume", "--cores", "3")
        assert count() == before
        assert proc.stdout.count("reused verified LASTZ checkpoint") == 4
        assert dump_db(cached / "probe.matches.sqlite") == dump_db(serial / "probe.matches.sqlite")
        # Output corruption must cause exactly one recomputation.
        target = sorted(cached.glob("*.lastz"))[0]
        target.write_text("corrupt\n")
        run("cached", "--resume", "--cores", "2")
        assert count() == before + 1
        before = count()
        # Same-size input mutation, with mtime restored, still invalidates.
        source = sorted(contigs.glob("*.fasta"))[0]
        stat = source.stat()
        data = source.read_bytes()
        start = data.index(b"\n") + 1
        changed = data[:start] + bytes([ord('C') if data[start] != ord('C') else ord('A')]) + data[start + 1:]
        source.write_bytes(changed)
        os.utime(source, ns=(stat.st_atime_ns, stat.st_mtime_ns))
        run("cached", "--resume", "--cores", "2")
        assert count() == before + 1
        before = count()
        run("cached", "--resume", "--min-identity", "81")
        assert count() == before + 4
        before = count()
        probes.write_text(probes.read_text() + "\n")
        run("cached", "--resume", "--min-identity", "81")
        assert count() == before + 4
        before = count()
        tool.write_text(tool.read_text() + "\n# changed executable\n")
        run("cached", "--resume", "--min-identity", "81")
        assert count() == before + 4
        # Missing completion marker cannot be inferred from a .lastz alone.
        next(cached.glob("*.complete")).unlink()
        before = count()
        run("cached", "--resume", "--min-identity", "81")
        assert count() == before + 1
        next(cached.glob("*.complete")).write_bytes(b"\xff")
        before = count()
        run("cached", "--resume", "--min-identity", "81")
        assert count() == before + 1
        # Failed samples do not publish partial files; completed samples survive.
        fail_name = sorted(contigs.glob("*.fasta"))[1].stem
        failed, _ = run("failed", "--resume", fail=True, extra_env={"TEST_LASTZ_FAIL": fail_name})
        assert len(list(failed.glob("*.complete"))) == 1
        assert not (failed / (fail_name + ".lastz")).exists()
        assert not list(failed.glob("*.partial"))
        before = count()
        run("failed", "--resume", "--cores", "2")
        assert count() == before + 3
        assert len(list(failed.glob("*.complete"))) == 4
        # Argument errors must precede destructive output handling.
        run("serial", "--force", "--cores", "0", fail=True)
        assert (serial / "probe.matches.sqlite").exists()
        run("cached", "--resume", "--force", fail=True)
        run("cached", "--resume", "--skip-alignment", fail=True)
    print("parallel/resume: ordered outputs, selective invalidation and failure recovery pass.")


if __name__ == "__main__":
    main()
