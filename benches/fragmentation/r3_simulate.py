#!/usr/bin/env python3
"""Simulate reads with a known V-shaped gradient for the R3 experiment.

Same model as benches/simulate.py (Pilea's Methods): replication initiates at
position 0, terminates at the midpoint, coverage decreases log2-linearly from
ori to ter, so the tent amplitude is exactly log2PTR. Reads are drawn from the
COMPLETE genome; they are then profiled against every reference (the scrambled
references preserve sequence content, so only the coordinate changes -- that is
the whole point of the fragmentation protocol).

One FASTQ per (log2PTR, depth) cell, fixed seeds:

    python3 r3_simulate.py ../genomes/Escherichia_coli_K12.fna -o $WORK/reads
"""
import argparse
import os
import sys

import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))
from simulate import load, sample_reads  # noqa: E402

LOG2PTRS = [0.0, 0.5, 1.0, 1.5]
DEPTHS = [1, 5, 10]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("fasta")
    ap.add_argument("-o", "--outdir", required=True)
    ap.add_argument("--seed", type=int, default=1000)
    a = ap.parse_args()

    os.makedirs(a.outdir, exist_ok=True)
    seq = load(a.fasta)
    for i, lp in enumerate(LOG2PTRS):
        for j, cov in enumerate(DEPTHS):
            out = os.path.join(a.outdir, f"lp{lp:g}_{cov:g}x.fq")
            if os.path.exists(out):
                continue
            rng = np.random.default_rng(a.seed + 100 * i + j)
            reads = sample_reads(seq, cov, lp, rng)
            with open(out, "w") as fh:
                for k, r in enumerate(reads):
                    fh.write(f"@r{k}\n{r}\n+\n{'I' * len(r)}\n")
            print(f"{out}: {len(reads):,} reads (log2PTR={lp:g}, {cov:g}x)")


if __name__ == "__main__":
    main()
