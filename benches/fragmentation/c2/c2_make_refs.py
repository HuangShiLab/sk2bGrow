#!/usr/bin/env python3
"""C2 Task 2.1/2.3: genuinely incomplete MAG references (not just cut ones).

fragment.py preserves every base; a real MAG does not. This builds the
incomplete + contaminated variants of the 100-contig K-12 draft:

  compNN            drop a random (100-NN)% of the frag100 contigs (seed 0)
  compNN_cMM_<src>  the compNN draft plus MM foreign contigs spliced in from
                    a draft of <src> (cut with the identical fragment.py
                    protocol, seed 0; contigs chosen with rng seed 0)

Self contigs keep their fragment.py names (contig_XXXX) so the layout truth
still applies to them; foreign contigs are named fcontig_<src>_XXXX.

Also a --frag16 mode for Task 2.3: every genome in the 16-genome set cut
into 100 contigs (longest record; other records, e.g. plasmids, appended
unchanged), written with the ORIGINAL file stems so genome names in the
index match the complete-genome simulation truth.

Seeds are fixed (0) everywhere so the numbers are reproducible; the
fragment.py protocol (lognormal mu=0 sigma=1, shuffle, p(flip)=0.5) is
reused verbatim via import.
"""
import argparse
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from fragment import read_fasta, cut_points  # noqa: E402

COMP = bytes.maketrans(b"ACGTNacgtn", b"TGCANtgcan")


def write_fasta(records, path):
    os.makedirs(os.path.dirname(path) or ".", exist_ok=True)
    with open(path, "w") as fh:
        for name, s in records:
            fh.write(f">{name}\n")
            for j in range(0, len(s), 70):
                fh.write(s[j:j + 70] + "\n")


def frag_records(seq, n, seed):
    """Cut `seq` exactly the way fragment.py main() does."""
    rng = np.random.default_rng(seed)
    lens = cut_points(len(seq), n, rng)
    starts = np.r_[0, np.cumsum(lens)[:-1]]
    order = rng.permutation(n)
    flip = rng.random(n) < 0.5
    out = []
    for out_i, i in enumerate(order):
        s = seq[starts[i]:starts[i] + lens[i]]
        if flip[i]:
            s = s.translate(COMP)[::-1]
        out.append((f"contig_{out_i:04d}", s))
    return out


def load_frag100(path):
    recs = list(read_fasta(path))
    return [(n, s) for n, s in recs]


def foreign_contigs(genomes_dir, src_stem, k, seed, n=100):
    """k contigs from a fragment.py-protocol draft of the longest record of
    <src_stem>.fna (same seed/protocol as the self draft)."""
    path = os.path.join(genomes_dir, src_stem + ".fna")
    recs = sorted(read_fasta(path), key=lambda r: -len(r[1]))
    draft = frag_records(recs[0][1], n, seed)
    rng = np.random.default_rng(seed)
    take = rng.choice(len(draft), size=k, replace=False)
    return [(f"fcontig_{src_stem}_{draft[i][0]}", draft[i][1]) for i in take]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--frag", help="frag100 FASTA (from fragment.py)")
    ap.add_argument("--genomes-dir", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--seed", type=int, default=0)
    ap.add_argument("-n", "--n-contigs", type=int, default=100)
    ap.add_argument("--frag16", action="store_true",
                    help="fragment every genome in --genomes-dir instead")
    a = ap.parse_args()

    if a.frag16:
        os.makedirs(a.out, exist_ok=True)
        for fn in sorted(os.listdir(a.genomes_dir)):
            if not fn.endswith(".fna"):
                continue
            recs = list(read_fasta(os.path.join(a.genomes_dir, fn)))
            longest = max(range(len(recs)), key=lambda i: len(recs[i][1]))
            draft = frag_records(recs[longest][1], a.n_contigs, a.seed)
            kept = [r for i, r in enumerate(recs) if i != longest]
            write_fasta(draft + kept, os.path.join(a.out, fn))
            print(f"frag16 {fn}: {len(draft)} contigs from longest record "
                  f"({len(recs[longest][1]):,} bp) + {len(kept)} records kept")
        return

    frag = load_frag100(a.frag)
    rng = np.random.default_rng(a.seed)
    manifest = ["condition\tn_contigs\tn_self\tn_dropped\tbp_self\tcompleteness"
                "\tforeign_src\tn_foreign\tforeign_bp"]
    total_bp = sum(len(s) for _, s in frag)

    # completeness drop sets are nested decisions drawn once from seed 0:
    # one permutation, prefixes dropped -> comp90 subset of comp75 etc. is NOT
    # wanted (independent draws are the fair test), so draw per condition.
    variants = []
    for comp in (90, 75, 50):
        n_drop = 100 - comp
        variants.append((f"comp{comp}", comp, None, 0))
    variants += [("comp75_c10_O157H7", 75, "Escherichia_coli_O157H7", 10),
                 ("comp50_c10_O157H7", 50, "Escherichia_coli_O157H7", 10),
                 ("comp75_c10_Salmonella", 75, "Salmonella_enterica_LT2", 10)]

    for name, comp, fsrc, nfor in variants:
        n_drop = 100 - comp
        drop = set(rng.choice(len(frag), size=n_drop, replace=False).tolist())
        kept = [r for i, r in enumerate(frag) if i not in drop]
        bp_self = sum(len(s) for _, s in kept)
        foreign = []
        if fsrc:
            foreign = foreign_contigs(a.genomes_dir, fsrc, nfor, a.seed)
        write_fasta(kept + foreign, os.path.join(a.out, name + ".fna"))
        fbp = sum(len(s) for _, s in foreign)
        manifest.append(
            f"{name}\t{len(kept) + len(foreign)}\t{len(kept)}\t{n_drop}"
            f"\t{bp_self}\t{bp_self / total_bp:.4f}"
            f"\t{fsrc or '-'}\t{len(foreign)}\t{fbp}")
        print(f"{name}: {len(kept)} self ({bp_self:,} bp, "
              f"{100 * bp_self / total_bp:.1f}%) + {len(foreign)} foreign")

    with open(os.path.join(a.out, "refs_manifest.tsv"), "w") as fh:
        fh.write("\n".join(manifest) + "\n")


if __name__ == "__main__":
    main()
