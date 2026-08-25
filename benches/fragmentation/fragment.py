#!/usr/bin/env python3
"""Split a complete genome into draft contigs, the way Pilea's Fig 3 does.

Contig lengths are lognormal (mu = 0, sigma = 1) scaled to the genome length,
the order is shuffled and each contig is independently reverse-complemented with
probability 0.5. That destroys exactly what the V-shape fit depends on -- the
genomic coordinate -- while leaving the sequence content untouched, so any
change in accuracy is attributable to fragmentation alone.

    python3 fragment.py ecoli.fna -n 100 --seed 0 -o ecoli_frag100.fna

Writes a `.layout.tsv` beside the FASTA recording where each contig came from,
so a scaffolding result can be scored against the truth rather than only against
the downstream PTR.
"""
import argparse
import os

import numpy as np

COMP = bytes.maketrans(b"ACGTNacgtn", b"TGCANtgcan")


def read_fasta(path):
    name, chunks = None, []
    with open(path) as fh:
        for line in fh:
            if line.startswith(">"):
                if name is not None:
                    yield name, "".join(chunks)
                name, chunks = line[1:].strip(), []
            else:
                chunks.append(line.strip())
    if name is not None:
        yield name, "".join(chunks)


def cut_points(total, n, rng):
    """n lognormal lengths summing to `total`, none shorter than 500 bp (the
    indexer's own floor -- a contig it would discard is not a fair test)."""
    for _ in range(1000):
        w = rng.lognormal(0.0, 1.0, n)
        lens = np.maximum(500, np.floor(w / w.sum() * total).astype(np.int64))
        if lens.sum() <= total:
            lens[-1] += total - lens.sum()
            return lens
    raise RuntimeError("could not fit lognormal contig lengths; try fewer contigs")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("fasta")
    ap.add_argument("-n", "--n-contigs", type=int, default=100)
    ap.add_argument("--seed", type=int, default=0)
    ap.add_argument("-o", "--out", required=True)
    a = ap.parse_args()

    rng = np.random.default_rng(a.seed)
    records = list(read_fasta(a.fasta))
    if len(records) != 1:
        raise SystemExit(f"{a.fasta} has {len(records)} contigs; expected a complete genome")
    _, seq = records[0]

    lens = cut_points(len(seq), a.n_contigs, rng)
    starts = np.r_[0, np.cumsum(lens)[:-1]]
    order = rng.permutation(a.n_contigs)
    flip = rng.random(a.n_contigs) < 0.5

    os.makedirs(os.path.dirname(a.out) or ".", exist_ok=True)
    layout = open(os.path.splitext(a.out)[0] + ".layout.tsv", "w")
    layout.write("contig\ttrue_start\tlength\treversed\n")
    with open(a.out, "w") as fh:
        for out_i, i in enumerate(order):
            s = seq[starts[i]:starts[i] + lens[i]]
            if flip[i]:
                s = s.translate(COMP)[::-1]
            fh.write(f">contig_{out_i:04d}\n")
            for j in range(0, len(s), 70):
                fh.write(s[j:j + 70] + "\n")
            layout.write(f"contig_{out_i:04d}\t{starts[i]}\t{lens[i]}\t{int(flip[i])}\n")
    layout.close()
    print(f"{a.out}: {a.n_contigs} contigs, {lens.min():,}-{lens.max():,} bp "
          f"(N50 {n50(lens):,}), {flip.sum()} reversed")


def n50(lens):
    s = np.sort(lens)[::-1]
    return int(s[np.searchsorted(np.cumsum(s), s.sum() / 2)])


if __name__ == "__main__":
    main()
