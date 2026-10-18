#!/usr/bin/env python3
"""Warm-cache stage benchmarks; exact output checks, GNU time, real MAFFT.
Pass MAFFT_BINARIES in the environment when using an extracted MAFFT package.
"""
import argparse
import hashlib
import json
import os
import platform
import random
import statistics
import subprocess
import tempfile
import time
from pathlib import Path


def hashes(directory):
    return {p.name: hashlib.sha256(p.read_bytes()).hexdigest()
            for p in sorted(directory.iterdir()) if p.is_file()}


def measure(command, root, env):
    start = time.perf_counter()
    proc = subprocess.run(['/usr/bin/time', '-f', '%U %S %M', '-o', str(root / 'usage'), *command],
                          capture_output=True, env=env)
    elapsed = time.perf_counter() - start
    if proc.returncode:
        raise RuntimeError(proc.stderr.decode())
    user, system, rss = map(float, (root / 'usage').read_text().split())
    return dict(seconds=elapsed, user_seconds=user, system_seconds=system, max_rss_kib=rss), proc.stdout


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--baseline', type=Path, required=True)
    parser.add_argument('--candidate', type=Path, required=True)
    parser.add_argument('--mafft', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--repeats', type=int, default=3)
    parser.add_argument('--concat-only', action='store_true')
    args = parser.parse_args()
    env = os.environ.copy()
    results = {}
    with tempfile.TemporaryDirectory(prefix='phyluce-alignment-benchmark-') as temp:
        root = Path(temp)
        conf = root / 'phyluce.conf'
        conf.write_text('[binaries]\nmafft:' + str(args.mafft.resolve()) + '\n')
        env['PHYLUCE_CONFIG'] = str(conf)
        inputs = root / 'inputs'
        inputs.mkdir()
        # 100 MB full matrix, 80% occupancy, changing taxon membership/order.
        sequence = 'ACGT' * 125
        for locus in range(1000):
            with (inputs / f'locus-{locus:04}.fasta').open('w') as out:
                for taxon in reversed(range(200)):
                    if (taxon + locus) % 5:
                        out.write(f'>taxon_{taxon}\n{sequence}\n')
        rng = random.Random(7941)
        for shape in ['uniform', 'skewed']:
            with (root / f'{shape}.fasta').open('w') as out:
                for locus in range(20):
                    length = 400 if shape == 'uniform' or locus < 16 else 1800
                    taxa = 8 if length == 400 else 24
                    template = ''.join(rng.choices('ACGT', k=length))
                    for taxon in range(taxa):
                        seq = ''.join(rng.choice('ACGT') if rng.random() < .06 else c for c in template)
                        out.write(f'>taxon_{taxon} |locus-{locus:03}\n{seq}\n')
        cases = []
        for fmt in ['nexus', 'phylip']:
            cases.append((f'concat-{fmt}', ['align', 'concatenate-alignments', '--alignments', str(inputs),
                                           '--input-format', 'fasta', '--' + fmt],
                          {'baseline': [], 'candidate-disk': ['--memory-limit-mb', '0'],
                           'candidate-100mib': ['--memory-limit-mb', '100'], 'candidate': []}))
        for shape in ['uniform', 'skewed']:
            for trim in [False, True]:
                cases.append((f'mafft-{shape}-' + ('trim' if trim else 'raw'),
                              ['align', 'seqcap-align', '--input', str(root / f'{shape}.fasta'),
                               '--taxa', '8', '--incomplete-matrix', '--cores', '4'] + ([] if trim else ['--no-trim']),
                              {'baseline': [], 'candidate': []}))
        if args.concat_only:
            cases = [case for case in cases if case[0].startswith('concat-')]
        for name, options, variants in cases:
            samples = {key: [] for key in variants}
            reference = None
            for repeat in range(args.repeats):
                order = list(variants)
                if repeat % 2:
                    order.reverse()
                for key in order:
                    run = root / f'{name}-{repeat}-{key}'
                    run.mkdir()
                    output = run / 'result'
                    binary = args.baseline if key == 'baseline' else args.candidate
                    measurement, stdout = measure([str(binary.resolve()), *options, '--output', str(output), *variants[key]], run, env)
                    fingerprint = (hashes(output), stdout)
                    if reference is None:
                        reference = fingerprint
                    assert fingerprint == reference, (name, repeat, key, 'output differs')
                    samples[key].append(measurement)
            results[name] = {'output_identical': True, 'samples': samples,
                             'median': {key: {metric: statistics.median(s[metric] for s in values)
                                              for metric in values[0]} for key, values in samples.items()}}
            print(name, results[name]['median'], flush=True)
        version = subprocess.run([str(args.mafft.resolve()), '--version'], env=env, capture_output=True, text=True)
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(json.dumps({'platform': platform.platform(), 'mafft_version': version.stderr.strip() or version.stdout.strip(),
                                          'repeats': args.repeats, 'warm_cache': True,
                                          'binary_sha256': {k: hashlib.sha256(getattr(args, k).read_bytes()).hexdigest() for k in ['baseline', 'candidate']},
                                          'results': results}, indent=2) + '\n')


if __name__ == '__main__':
    main()
