#!/usr/bin/env python3
"""Reproducible synthetic stage benchmarks (Linux GNU time; no external aligner).
Outputs JSON measurements and verifies equivalent FASTA/SQLite/CSV outputs.
No cache flushing: numbers describe repeated, predominantly warm-cache runs.
"""
import argparse
import hashlib
import json
import platform
import shutil
import sqlite3
import statistics
import subprocess
import tempfile
import time
from pathlib import Path


def digest(path):
    with path.open("rb") as file:
        return hashlib.file_digest(file, "sha256").hexdigest()


def measure(command, directory):
    stats = directory / "usage.txt"
    start = time.perf_counter()
    with (directory / "stdout.txt").open("w") as stdout:
        proc = subprocess.run(["/usr/bin/time", "-f", "%U %S %M %I %O", "-o", str(stats), *command],
                              stdout=stdout, stderr=subprocess.PIPE, text=True)
    elapsed = time.perf_counter() - start
    if proc.returncode:
        raise RuntimeError(proc.stderr)
    user, system, rss, reads, writes = map(float, stats.read_text().split())
    return dict(wall_seconds=elapsed, user_seconds=user, system_seconds=system,
                max_rss_kib=int(rss), filesystem_inputs=int(reads), filesystem_outputs=int(writes))


def db_digest(path):
    with sqlite3.connect(path) as conn:
        rows = [(table, conn.execute(f'SELECT * FROM "{table}" ORDER BY uce').fetchall())
                for table in ["matches", "match_map"]]
    return hashlib.sha256(repr(rows).encode()).hexdigest()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--baseline", type=Path, required=True)
    parser.add_argument("--candidate", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--records", type=int, default=250000)
    parser.add_argument("--length", type=int, default=1000)
    parser.add_argument("--samples", type=int, default=24)
    parser.add_argument("--loci", type=int, default=3000)
    parser.add_argument("--cores", type=int, default=4)
    parser.add_argument("--repeats", type=int, default=3)
    args = parser.parse_args()
    assert min(args.records, args.length, args.samples, args.loci, args.cores, args.repeats) > 0
    binaries = {k: str(getattr(args, k).resolve()) for k in ["baseline", "candidate"]}
    results = {}
    with tempfile.TemporaryDirectory(prefix="phyluce-performance-") as temp:
        root = Path(temp)
        contigs = root / "extract-contigs"
        contigs.mkdir()
        sequence = ("ACGT" * ((args.length + 3) // 4))[:args.length]
        fasta = contigs / "sample.contigs.fasta"
        with fasta.open("w") as out:
            for i in range(args.records):
                out.write(f">NODE_{i}_length_{args.length}_cov_1.0\n{sequence}\n")
        indices = sorted(set([*range(min(100, args.records)), *range(max(0, args.records - 100), args.records)]))
        config = root / "taxa.conf"
        config.write_text("[Organisms]\nsample\n[Loci]\n" + "".join(f"uce-{i}\n" for i in indices))
        database = root / "extract.sqlite"
        with sqlite3.connect(database) as conn:
            conn.execute("CREATE TABLE match_map (uce TEXT PRIMARY KEY, sample TEXT)")
            conn.executemany("INSERT INTO match_map VALUES (?, ?)",
                [(f"uce-{i}", f"NODE_{i}(+)") for i in indices])

        match_contigs = root / "match-contigs"
        match_contigs.mkdir()
        probes = root / "probes.fasta"
        with probes.open("w") as out:
            for i in range(args.loci):
                out.write(f">uce-{i}_p1\n{'A' * 120}\n")
        source = root / "source.lastz"
        with source.open("w") as out:
            for i in range(args.loci):
                out.write(f"12000\tNODE_{i}_length_120_cov_1.0\t+\t0\t120\t120\tuce-{i}_p1\t+\t0\t120\t120\t{'.' * 120}\t120M\t120/120\t100.0%\t120/120\t100.0%\n")
        # Unique locus/contig hits, with a duplicate-probe row to exercise deduplication.
        for i in range(args.samples):
            with (match_contigs / f"sample{i:04d}.fasta").open("w") as out:
                for j in range(args.loci):
                    out.write(f">NODE_{j}_length_120_cov_1.0\n{'A' * 120}\n")
        stages = {"extraction": [], "match_result_processing": []}
        reference = {}
        for repetition in range(args.repeats):
            # Alternate order to reduce a systematic cache advantage.
            variants = ["baseline", "candidate"] if repetition % 2 == 0 else ["candidate", "baseline"]
            for variant in variants:
                run_dir = root / f"{variant}-{repetition}"
                run_dir.mkdir()
                output = run_dir / "extracted.fasta"
                command = [binaries[variant], "assembly", "get-fastas-from-match-counts",
                           "--contigs", str(contigs), "--locus-db", str(database),
                           "--match-count-output", str(config), "--output", str(output)]
                sample = measure(command, run_dir) | {"variant": variant}
                key = digest(output)
                assert key == reference.setdefault("extract", key)
                stages["extraction"].append(sample)
                for cores in ([1] if variant == "baseline" else sorted(set([1, args.cores]))):
                    out = run_dir / f"match-{cores}"
                    out.mkdir()
                    for i in range(args.samples):
                        shutil.copyfile(source, out / f"sample{i:04d}.lastz")
                    csv = run_dir / f"match-{cores}.csv"
                    command = [binaries[variant], "assembly", "match-contigs-to-probes",
                               "--contigs", str(match_contigs), "--probes", str(probes),
                               "--output", str(out), "--csv", str(csv), "--skip-alignment"]
                    if variant == "candidate": command += ["--cores", str(cores)]
                    sample = measure(command, run_dir) | {"variant": variant, "cores": cores}
                    key = (db_digest(out / "probe.matches.sqlite"), digest(csv))
                    assert key == reference.setdefault("match", key)
                    stages["match_result_processing"].append(sample)
        for stage, samples in stages.items():
            groups = {}
            for sample in samples:
                group = sample["variant"] + (f'-cores-{sample["cores"]}' if "cores" in sample else "")
                groups.setdefault(group, []).append(sample)
            results[stage] = {group: {
                "median_wall_seconds": statistics.median(s["wall_seconds"] for s in samples),
                "median_max_rss_kib": statistics.median(s["max_rss_kib"] for s in samples),
                "runs": samples,
            } for group, samples in groups.items()}
        results["metadata"] = {"platform": platform.platform(), "settings": {k: str(v) if isinstance(v, Path) else v for k, v in vars(args).items()},
                               "extraction_input_bytes": fasta.stat().st_size,
                               "binary_sha256": {k: digest(Path(v)) for k, v in binaries.items()},
                               "notes": "Synthetic warm-cache stage measurements; LASTZ is NOT executed. Outputs verified equivalent."}
    args.output.write_text(json.dumps(results, indent=2) + "\n")
    print(json.dumps({stage: {group: {k: v for k, v in stats.items() if k != "runs"} for group, stats in groups.items()}
                      for stage, groups in results.items() if stage != "metadata"}, indent=2))


if __name__ == "__main__":
    main()
