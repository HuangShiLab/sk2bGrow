#!/usr/bin/env python3
"""Tent-gradient read simulator with analytically known per-anchor truth.

Experiment C7 (RESEARCH_PLAN.md Part 6): reads are drawn from a known V-shaped
coverage gradient, so the *expected* count of every anchor is known exactly —
not just the genome-level PTR. That is what lets the window-rate layer
(ZTP/ZTNB, ``python/sk2bgrow/ztp.py``) be checked against truth instead of
against a downstream fit.

Model: log2 coverage falls linearly from the origin to the terminus,
``log2 c(x) = const - log2_ptr * d(x) / (L/2)`` with ``d`` the circular
distance. Read starts ``s`` are drawn with probability proportional to
``2 ** (-log2_ptr * d(s + READ_LEN/2) / (L/2))`` — because the gradient is
linear in log space, evaluating the weight at the read midpoint makes the
resulting per-base coverage tent exact (the average of a linear function over
the read span is its value at the centre).

A read is counted against an anchor iff it covers the anchor's whole tag,
``contig[position : position + tag_len]`` (digest.rs convention, forward
strand), i.e. start ``s`` in ``[position + tag_len - READ_LEN, position]``.
The expected count is therefore exactly::

    lambda(anchor) = n_reads * sum_{s=lo}^{hi} p(s),   p = w / sum(w)

which :func:`anchor_expected` computes with a prefix sum. ``analyze.py``
rebuilds ``w`` with the same functions from the recorded parameters (plus the
per-sample efficiency vector for the overdispersed arm), so the truth and the
reads can never disagree about the model.

Overdispersed arm (``--sigma-eff > 0``): each 1 kb bin gets a fixed lognormal
efficiency multiplier (unit mean), so within-window counts are genuinely
overdispersed and the ZTNB branch becomes identifiable. The bin vector is
saved to ``--eff-out``; truth reconstruction multiplies the same weights by it.
"""
from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np

READ_LEN = 150
#: MG1655 oriC (DoriC), the same constant as ``sk2bgrow.simulate.ECOLI_ORI``.
ORI = 3_923_883
#: Efficiency bins are 1 kb.
EFF_BIN = 1000


def load_fasta(path: str | Path) -> bytes:
    """Longest record, upper-case — same convention as benches/simulate.py."""
    best, cur = b"", []
    with open(path, "rb") as fh:
        for line in fh:
            if line.startswith(b">"):
                if len(b"".join(cur)) > len(best):
                    best = b"".join(cur)
                cur = []
            else:
                cur.append(line.strip())
    if len(b"".join(cur)) > len(best):
        best = b"".join(cur)
    return best.upper()


def start_weights(genome_len: int, ori: float, log2_ptr: float, eff: np.ndarray | None = None) -> np.ndarray:
    """Unnormalised sampling weight of every valid read start.

    ``eff`` (optional) is one multiplier per ``EFF_BIN`` bp bin, applied by the
    bin of the read start. Deterministic given the arguments, so the analysis
    layer reconstructs the identical array.
    """
    starts = np.arange(genome_len - READ_LEN + 1, dtype=np.float64)
    mid = starts + READ_LEN / 2.0
    d = np.mod(mid - float(ori), float(genome_len))
    d = np.minimum(d, genome_len - d)
    w = np.exp2(-float(log2_ptr) * d / (genome_len / 2.0))
    if eff is not None:
        bins = np.minimum(starts.astype(np.int64) // EFF_BIN, len(eff) - 1)
        w = w * eff[bins]
    return w


def anchor_expected(w: np.ndarray, n_reads: int, positions: np.ndarray, tag_lens: np.ndarray) -> np.ndarray:
    """Exact expected count per anchor under start-weight vector ``w``.

    ``positions``/``tag_lens`` are per anchor (0-based forward-strand tag start
    and tag length). A read covers the tag iff its start lies in
    ``[position + tag_len - READ_LEN, position]`` intersected with the valid
    start range, so anchors nearer than one read length to the linear contig
    ends get their (smaller) true expectation rather than a special case.
    """
    genome_len = len(w) + READ_LEN - 1
    cdf = np.concatenate(([0.0], np.cumsum(w)))
    pos = np.asarray(positions, dtype=np.int64)
    tl = np.asarray(tag_lens, dtype=np.int64)
    lo = np.maximum(0, pos + tl - READ_LEN)
    hi = np.minimum(pos, genome_len - READ_LEN)
    lam = n_reads * (cdf[hi + 1] - cdf[lo]) / cdf[-1]
    return np.where(lo > hi, 0.0, lam)


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--genome", required=True)
    ap.add_argument("--depth", type=float, required=True, help="mean per-base coverage")
    ap.add_argument("--log2ptr", type=float, required=True, help="tent amplitude = true log2(PTR)")
    ap.add_argument("--seed", type=int, required=True)
    ap.add_argument("--sigma-eff", type=float, default=0.0, help="lognormal sd of per-kb efficiency (0 = Poisson arm)")
    ap.add_argument("--ori", type=int, default=ORI)
    ap.add_argument("--out", required=True, help="output FASTQ")
    ap.add_argument("--eff-out", help="where to save the efficiency vector (npy); required with --sigma-eff > 0")
    args = ap.parse_args()

    seq = load_fasta(args.genome)
    L = len(seq)
    rng = np.random.default_rng(args.seed)

    eff = None
    if args.sigma_eff > 0:
        if not args.eff_out:
            ap.error("--sigma-eff requires --eff-out")
        eff = rng.lognormal(-0.5 * args.sigma_eff**2, args.sigma_eff, L // EFF_BIN + 1)
        Path(args.eff_out).parent.mkdir(parents=True, exist_ok=True)
        np.save(args.eff_out, eff)

    w = start_weights(L, args.ori, args.log2ptr, eff)
    n_reads = int(round(args.depth * L / READ_LEN))
    starts = rng.choice(w.size, size=n_reads, replace=True, p=w / w.sum())
    starts.sort()

    Path(args.out).parent.mkdir(parents=True, exist_ok=True)
    with open(args.out, "w") as fh:
        for i, s in enumerate(starts):
            read = seq[s : s + READ_LEN].decode("ascii")
            fh.write(f"@r{i}\n{read}\n+\n{'I' * len(read)}\n")
    print(f"{args.out}: {n_reads:,} reads, depth {args.depth}x, log2PTR {args.log2ptr}, seed {args.seed}")


if __name__ == "__main__":
    main()
