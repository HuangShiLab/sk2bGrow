#!/usr/bin/env python3
"""References for the R3 QC-statistic experiment.

Two families over the same cut points (fragment.py's lognormal protocol,
seed 0):

  scr<N>    N contigs, shuffled and flipped -- the MAG case; identical to what
            fragment.py writes, regenerated here so both families share cuts.
  ord<N>    N contigs at the SAME cut points but in true order and orientation
            -- the "correctly ordered multi-contig scaffold" control. The QC
            statistic must stay quiet on these: fragmentation alone is not a
            failure, only a *scrambled* coordinate is.

The complete genome itself is the N=1 reference and is used as-is.

    python3 r3_refs.py ../genomes/Escherichia_coli_K12.fna -o $WORK/refs
"""
import argparse
import os

import numpy as np

from fragment import COMP, cut_points, n50, read_fasta


def write_fasta(path, chunks):
    with open(path, "w") as fh:
        for name, s in chunks:
            fh.write(f">{name}\n")
            for j in range(0, len(s), 70):
                fh.write(s[j:j + 70] + "\n")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("fasta")
    ap.add_argument("-o", "--outdir", required=True)
    ap.add_argument("--counts", default="2 5 10 20 50 100")
    ap.add_argument("--ordered", default="10 100",
                    help="which cut sets also get a correctly-ordered control")
    ap.add_argument("--seed", type=int, default=0)
    ap.add_argument("--rotate", type=float, default=0.0,
                    help="rotate the (circular) chromosome by this fraction of its "
                         "length first; the simulator puts ori at position 0, which "
                         "is always a fragment cut point -- rotating puts ori at a "
                         "generic interior position instead")
    ap.add_argument("--suffix", default="")
    a = ap.parse_args()

    os.makedirs(a.outdir, exist_ok=True)
    records = list(read_fasta(a.fasta))
    if len(records) != 1:
        raise SystemExit(f"{a.fasta} has {len(records)} contigs; expected a complete genome")
    _, seq = records[0]
    if a.rotate:
        k = int(len(seq) * a.rotate)
        seq = seq[k:] + seq[:k]
        complete = os.path.join(a.outdir, f"complete{a.suffix}.fna")
        write_fasta(complete, [(f"complete{a.suffix}", seq)])
        print(f"{complete}: rotated by {k:,} bp")

    counts = [int(x) for x in a.counts.split()]
    ordered = {int(x) for x in a.ordered.split()}
    for n in counts:
        # same draws as fragment.py: lengths, then permutation, then flips
        rng = np.random.default_rng(a.seed)
        lens = cut_points(len(seq), n, rng)
        starts = np.r_[0, np.cumsum(lens)[:-1]]
        order = rng.permutation(n)
        flip = rng.random(n) < 0.5

        scr = []
        for out_i, i in enumerate(order):
            s = seq[starts[i]:starts[i] + lens[i]]
            if flip[i]:
                s = s.translate(COMP)[::-1]
            scr.append((f"contig_{out_i:04d}", s))
        path = os.path.join(a.outdir, f"scr{n}{a.suffix}.fna")
        write_fasta(path, scr)
        print(f"{path}: {n} contigs shuffled/flipped (N50 {n50(lens):,})")

        if n in ordered:
            ord_chunks = [(f"contig_{i:04d}", seq[starts[i]:starts[i] + lens[i]])
                          for i in range(n)]
            path = os.path.join(a.outdir, f"ord{n}{a.suffix}.fna")
            write_fasta(path, ord_chunks)
            print(f"{path}: {n} contigs in true order/orientation")


if __name__ == "__main__":
    main()
