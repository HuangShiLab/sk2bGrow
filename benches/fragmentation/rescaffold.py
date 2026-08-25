#!/usr/bin/env python3
"""Score a `sk2bgrow scaffold` result and rebuild a coordinate-bearing FASTA.

Two jobs, because they answer different questions:

  --score   how close the inferred placement is to the truth in `.layout.tsv`
            (placement rate, orientation accuracy, coordinate error). This is
            scaffolding accuracy on its own terms.

  -o FASTA  re-emit the draft as a single pseudo-contig with every contig
            written at its inferred reference start, padded with N. The indexer
            rejects tags containing N, so the padding contributes no anchors and
            costs nothing except a slightly longer sequence; what it buys is that
            `index` then assigns global coordinates matching the scaffold, which
            is what the V-shape fit needs. Contigs the scaffold could not place
            are appended past the end rather than dropped.

    python3 rescaffold.py draft.fna scaf.scaffold.json --score -o scaffolded.fna
"""
import argparse
import json

import numpy as np
from scipy import stats

COMP = bytes.maketrans(b"ACGTNacgtn", b"TGCANtgcan")


def read_fasta(path):
    name, chunks = None, []
    with open(path) as fh:
        for line in fh:
            if line.startswith(">"):
                if name is not None:
                    yield name, "".join(chunks)
                name, chunks = line[1:].split()[0], []
            else:
                chunks.append(line.strip())
    if name is not None:
        yield name, "".join(chunks)


def circ_err(a, b, span):
    """Placement error on a circular chromosome: 3 kb from the origin is 3 kb,
    not `span - 3000`."""
    d = np.abs(a - b) % span
    return np.minimum(d, span - d)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("draft")
    ap.add_argument("scaffold_json")
    ap.add_argument("--layout", default=None, help="truth from fragment.py")
    ap.add_argument("--score", action="store_true")
    ap.add_argument("-o", "--out", default=None)
    ap.add_argument("--label", default="scaffolded")
    a = ap.parse_args()

    seqs = dict(read_fasta(a.draft))
    names = list(seqs)
    res = json.load(open(a.scaffold_json))
    placed = {p["contig_id"]: p for p in res["placements"]}

    if a.score:
        layout = a.layout or a.draft.rsplit(".", 1)[0] + ".layout.tsv"
        truth = {}
        with open(layout) as fh:
            next(fh)
            for line in fh:
                n, start, length, rev = line.split()
                truth[n] = (int(start), int(length), bool(int(rev)))
        span = sum(v[1] for v in truth.values())
        errs, ori_ok, bp_placed = [], 0, 0
        for cid, p in placed.items():
            name = names[cid]
            t_start, t_len, t_rev = truth[name]
            errs.append(circ_err(np.int64(p["ref_start"]), np.int64(t_start), span))
            ori_ok += int((p["orientation"] == "Reverse") == t_rev)
            bp_placed += t_len
        errs = np.array(errs, dtype=float)
        print(f"placed        {len(placed)}/{len(names)} contigs "
              f"({100 * bp_placed / span:.1f}% of bp)")
        print(f"orientation   {ori_ok}/{len(placed)} correct "
              f"({100 * ori_ok / max(1, len(placed)):.1f}%)")
        print(f"start error   median {np.median(errs):,.0f} bp   "
              f"90th pct {np.percentile(errs, 90):,.0f} bp   max {errs.max():,.0f} bp")

        # Raw error over-states the damage. The estimator searches for the
        # origin rather than assuming it, so a rigid rotation of the whole
        # scaffold is invisible to the fit; scaffolding against a relative
        # whose origin sits elsewhere produces exactly that. Report the error
        # after removing the best circular shift, and the rank correlation of
        # inferred against true order, which is what the V-fit actually needs.
        inferred = np.array([placed[c]["ref_start"] for c in placed], dtype=float)
        true = np.array([truth[names[c]][0] for c in placed], dtype=float)

        def shifted(sh):
            d = np.abs((inferred - sh) % span - true) % span
            return np.minimum(d, span - d)

        best = min(np.linspace(0, span, 4001), key=lambda sh: np.median(shifted(sh)))
        e = shifted(best)
        rho = stats.spearmanr(inferred, true).statistic
        print(f"  after removing a {best:,.0f} bp rotation: median "
              f"{np.median(e):,.0f} bp   within 100 kb: {100 * (e < 100_000).mean():.0f}%")
        print(f"  contig order (Spearman, inferred vs true): {rho:.4f}")

    if a.out:
        span = max((p["ref_start"] + p["contig_len"] for p in placed.values()), default=0)
        canvas = bytearray(b"N" * int(span))
        for cid, p in placed.items():
            s = seqs[names[cid]].encode()
            if p["orientation"] == "Reverse":
                s = s.translate(COMP)[::-1]
            start = int(p["ref_start"])
            canvas[start:start + len(s)] = s[:len(canvas) - start]
        tail = b"".join(seqs[names[c]].encode() for c in res["unplaced"])
        seq = bytes(canvas) + tail
        with open(a.out, "w") as fh:
            fh.write(f">{a.label}\n")
            for j in range(0, len(seq), 70):
                fh.write(seq[j:j + 70].decode() + "\n")
        print(f"wrote {a.out}: {len(seq):,} bp "
              f"({100 * seq.count(b'N') / len(seq):.1f}% N)")


if __name__ == "__main__":
    main()
