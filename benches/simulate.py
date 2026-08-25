#!/usr/bin/env python3
"""Multi-strain community simulator, following Pilea's Methods verbatim.

    "We assumed that replication initiated at position zero (origin) and
     terminated at the midpoint (terminus) of the assembly. Coverage was
     simulated to decrease log2-linearly from the origin to the terminus by
     adjusting the sampling weights accordingly (V-shaped coverage profile).
     The values of log2(PTR) were drawn randomly from a uniform [0,2]."

Reads are written single-end; Pilea concatenated R1/R2 for the tools that do not
accept pairs anyway, so this removes a variable rather than adding one.
"""
import argparse, gzip, os
from pathlib import Path
import numpy as np

BIN = 1000          # position bins for the sampling weights
READ_LEN = 150


def load(fa):
    """Longest record only — the replication profile is defined on the chromosome."""
    best, cur = b'', []
    with open(fa, 'rb') as fh:
        for line in fh:
            if line.startswith(b'>'):
                s = b''.join(cur)
                if len(s) > len(best):
                    best = s
                cur = []
            else:
                cur.append(line.strip())
    s = b''.join(cur)
    return (s if len(s) > len(best) else best).upper()


def sample_reads(seq, cov, log2ptr, rng):
    """Draw reads with a V-shaped coverage profile around ori = 0."""
    L = len(seq)
    n_reads = int(cov * L / READ_LEN)
    nb = max(L // BIN, 2)
    centres = (np.arange(nb) + 0.5) * BIN
    d = np.minimum(centres, L - centres) / (L / 2.0)     # 0 at ori, 1 at ter
    w = np.exp2(-log2ptr * d)
    w /= w.sum()
    counts = rng.multinomial(n_reads, w)
    out = []
    rc = str.maketrans('ACGTN', 'TGCAN')
    for b in np.nonzero(counts)[0]:
        k = counts[b]
        starts = rng.integers(b * BIN, min((b + 1) * BIN, max(L - READ_LEN, 1)), size=k)
        for s in starts:
            r = seq[s:s + READ_LEN].decode('ascii', 'ignore')
            if len(r) < READ_LEN:
                continue
            if rng.random() < 0.5:
                r = r.translate(rc)[::-1]
            out.append(r)
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--genomes', default='genomes')
    ap.add_argument('--n-strains', type=int, required=True)
    ap.add_argument('--coverage', type=float, required=True)
    ap.add_argument('--seed', type=int, default=0)
    ap.add_argument('--out', required=True)
    ap.add_argument('--truth', required=True)
    a = ap.parse_args()

    rng = np.random.default_rng(a.seed)
    fas = sorted(Path(a.genomes).glob('*.fna'))
    pick = rng.choice(len(fas), size=min(a.n_strains, len(fas)), replace=False)

    rows, n = [], 0
    with open(a.out, 'w') as fh:
        for i in pick:
            fa = fas[i]
            seq = load(fa)
            log2ptr = float(rng.uniform(0, 2))
            for r in sample_reads(seq, a.coverage, log2ptr, rng):
                fh.write(f'@r{n}\n{r}\n+\n{"I" * len(r)}\n')
                n += 1
            rows.append((fa.stem, len(seq), log2ptr))
    with open(a.truth, 'w') as fh:
        fh.write('genome\tlength\ttrue_log2ptr\n')
        for g, L, p in rows:
            fh.write(f'{g}\t{L}\t{p:.6f}\n')
    print(f'{a.out}: {n:,} reads, {len(rows)} strains @ {a.coverage}x')


if __name__ == '__main__':
    main()
