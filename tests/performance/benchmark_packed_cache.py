#!/usr/bin/env python3
"""Compare packed concatenation cache with a supplied baseline; GNU time required."""
import argparse
import hashlib
import json
import os
import platform
import statistics
import tempfile
import threading
from pathlib import Path
from benchmark_alignment import hashes, measure


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    for option in ['baseline', 'candidate', 'output']:
        parser.add_argument('--' + option, required=True, type=Path)
    parser.add_argument('--repeats', type=int, default=5)
    parser.add_argument('--input-format', choices=['fasta', 'nexus'], default='fasta')
    parser.add_argument('--cases', nargs='+', choices=['standard', 'short-loci', 'sparse'], default=['standard', 'short-loci', 'sparse'])
    parser.add_argument('--track-temp-disk', action='store_true')
    args = parser.parse_args()
    if args.repeats < 1:
        parser.error('--repeats must be positive')
    results = {}
    with tempfile.TemporaryDirectory(prefix='phyluce-packed-') as temp:
        root = Path(temp)
        for case, loci, taxa, length, occupancy in [('standard', 1000, 200, 500, 8),
                                                  ('short-loci', 5000, 200, 40, 8),
                                                  ('sparse', 1000, 200, 500, 1)]:
            if case not in args.cases:
                continue
            inputs = root / case
            inputs.mkdir()
            for locus in range(loci):
                with (inputs / f'locus-{locus:05}.{args.input_format}').open('w') as out:
                    if args.input_format == 'nexus':
                        out.write(f'#NEXUS\nbegin data;\ndimensions ntax={taxa * occupancy // 10} nchar={length};\nformat datatype=dna missing=? gap=-;\nmatrix\n')
                    for taxon in reversed(range(taxa)):
                        if (taxon + locus) % 10 < occupancy:
                            sequence = ''.join('ACGTN-?'[(taxon + locus + i) % 7] for i in range(length))
                            if args.input_format == 'fasta':
                                out.write(f'>taxon_{taxon:03}_long_unique_name\n{sequence}\n')
                            else:
                                out.write(f'taxon_{taxon:03}_long_unique_name {sequence}\n')
                    if args.input_format == 'nexus':
                        out.write(';\nend;\n')
            for fmt in ['nexus', 'phylip']:
                # Exercise default, forced disk, and a budget near compact-cache size.
                variants = {'baseline': (args.baseline, 256), 'candidate': (args.candidate, 256),
                            'candidate-disk': (args.candidate, 0), 'baseline-disk': (args.baseline, 0),
                            'baseline-48mib': (args.baseline, 48), 'candidate-48mib': (args.candidate, 48),
                            'baseline-budget': (args.baseline, 32), 'candidate-budget': (args.candidate, 32)}
                samples = {key: [] for key in variants}
                reference = None
                for repeat in range(args.repeats):
                    order = list(variants)
                    if repeat % 2:
                        order.reverse()
                    for key in order:
                        binary, budget = variants[key]
                        run = root / f'{case}-{fmt}-{repeat}-{key}'
                        run.mkdir()
                        output = run / 'result'
                        scratch = run / 'scratch'
                        scratch.mkdir()
                        env = os.environ.copy()
                        env['TMPDIR'] = str(scratch)
                        stop = threading.Event()
                        peak = [0, 0]
                        def sample():
                            while not stop.is_set():
                                logical = allocated = 0
                                for path in scratch.iterdir():
                                    try:
                                        stat = path.stat()
                                    except FileNotFoundError:
                                        continue
                                    logical += stat.st_size
                                    allocated += stat.st_blocks * 512
                                peak[0] = max(peak[0], logical)
                                peak[1] = max(peak[1], allocated)
                                stop.wait(0.001)
                        watcher = threading.Thread(target=sample) if args.track_temp_disk else None
                        if watcher:
                            watcher.start()
                        try:
                            stats, stdout = measure([str(binary.resolve()), 'align', 'concatenate-alignments',
                                                 '--alignments', str(inputs), '--input-format', args.input_format,
                                                 '--output', str(output), '--' + fmt,
                                                 '--memory-limit-mb', str(budget)], run, env)
                        finally:
                            stop.set()
                            if watcher:
                                watcher.join()
                        assert not list(scratch.iterdir()), 'temporary files leaked'
                        if watcher:
                            stats.update(temp_peak_bytes=peak[0], temp_peak_allocated_bytes=peak[1])
                        fingerprint = (hashes(output), stdout)
                        if reference is None:
                            reference = fingerprint
                        assert fingerprint == reference, (case, fmt, repeat, key)
                        samples[key].append(stats)
                        # Avoid retaining gigabytes of repeated output across cases.
                        for file in output.iterdir():
                            file.unlink()
                median = {key: {metric: statistics.median(s[metric] for s in values)
                                for metric in values[0]} for key, values in samples.items()}
                results[f'{case}-{fmt}'] = dict(loci=loci, taxa=taxa, length=length,
                                                occupancy=occupancy / 10, output_identical=True,
                                                samples=samples, median=median)
                print(case, fmt, median, flush=True)
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(json.dumps(dict(platform=platform.platform(), repeats=args.repeats,
                                               warm_cache=True, temp_disk_sample_ms=1 if args.track_temp_disk else None, input_format=args.input_format, results=results,
                                               binary_sha256={key: hashlib.sha256(getattr(args, key).read_bytes()).hexdigest()
                                                              for key in ['baseline', 'candidate']}), indent=2) + '\n')


if __name__ == '__main__':
    main()
