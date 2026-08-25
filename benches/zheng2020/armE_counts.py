"""Arm E: rewrite Pilea's FracMinHash sketch as a sk2bGrow count table.

This is the missing cell of the 2 x 2 (sketch x estimator). Rather than
re-implementing our estimator on Pilea's data structures, it rewrites Pilea's
sketch into the count-table format the sk2bGrow statistics layer already
consumes, so the *entire* estimator -- GC correction, ZTP/ZTNB window rates with
standard errors, adaptive windows, shared-origin search, V-shape fit, fusion QC
-- runs on both sketches unchanged. The only difference between this arm and
sk2bGrow proper is which loci are counted.

Feeding raw window rates into the V-fit instead would not be a controlled
comparison: it returns log2PTR 2.65 against a measured 1.73, because what was
removed is the surrounding machinery rather than the fitting geometry.

Sketch positions are not stored in Pilea's .pdb, so they are recovered by
replaying `hash64` over the reference exactly as its `sketch.py` does.

Run with Pilea's interpreter (it imports pilea):

    pilea_env/bin/python armE_counts.py -r ecoli.fna -o counts sub/*.fq
"""
import argparse
import os
import pickle
import sys
from collections import defaultdict

import numpy as np
from pilea.io import parse_fastx_file
from pilea.kmc import hash64, count64
from pilea.sketch import GC_LUT, scan

FLAG_UNIQUE = (1 << 0) | (1 << 1)   # unique in genome | unique across db
FLAG_MULTICOPY = 1 << 2
HEADER = ("sample\tgenome_id\tgenome\tcontig_id\tposition\tglobal_position\tenzyme"
          "\tstrand\tflags\tlocal_gc\twindow_id\tcount")


def sketch_loci(fasta, k, maxhash):
    """{key: [(contig, pos, local_gc), ...]} plus per-contig lengths, computed
    the way pilea.sketch computes them."""
    loci, lens = defaultdict(list), []
    for cid, seq in enumerate(parse_fastx_file(fasta)):
        arr = np.frombuffer(seq, dtype=np.uint8)
        gc_prefix = np.r_[np.uint32(0), np.cumsum(GC_LUT[arr], dtype=np.uint32)]
        lens.append(len(seq))
        for i, key in hash64(seq, k, maxhash):
            loci[key].append((cid, i, scan(i + k // 2, gc_prefix)))
    return dict(loci), lens


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("reads", nargs="+")
    ap.add_argument("-r", "--reference", required=True)
    ap.add_argument("-g", "--genome", default="ecoli")
    ap.add_argument("-k", type=int, default=31)
    ap.add_argument("-s", type=int, default=250)
    ap.add_argument("-o", "--outdir", required=True)
    ap.add_argument("--cache", default=None, help="pickle of the reference sketch")
    a = ap.parse_args()

    maxhash = ((1 << 64) - 1) // a.s
    cache = a.cache or os.path.join(a.outdir, "sketch_loci.pkl")
    os.makedirs(a.outdir, exist_ok=True)
    if os.path.exists(cache):
        loci, lens = pickle.load(open(cache, "rb"))
    else:
        loci, lens = sketch_loci(a.reference, a.k, maxhash)
        pickle.dump((loci, lens), open(cache, "wb"))
    print(f"sketch: {len(loci):,} keys over {sum(lens):,} bp", file=sys.stderr)

    offsets, run = [], 0
    for length in lens:
        offsets.append(run)
        run += length

    for fq in a.reads:
        sample = os.path.basename(fq).rsplit(".", 1)[0]
        out = os.path.join(a.outdir, f"{sample}.counts.tsv")
        if os.path.exists(out):
            continue
        # count64 dedupes within a read (it builds a set per record). That is
        # Pilea's counting convention and must be preserved, so call it rather
        # than re-deriving counts from hash64.
        kmc = {}
        for record in parse_fastx_file(fq):
            count64(record, None, a.k, maxhash, kmc)
        with open(out, "w") as fh:
            fh.write(HEADER + "\n")
            for key, places in loci.items():
                # A k-mer occurring more than once in the genome has no single
                # locus; flag it multi-copy so the coverage model excludes it,
                # which is what Pilea's own unique/duplicate split does.
                flags = FLAG_MULTICOPY if len(places) > 1 else FLAG_UNIQUE
                cnt = int(kmc.get(key, 0))
                for cid, pos, gc in places:
                    fh.write(f"{sample}\t0\t{a.genome}\t{cid}\t{pos}"
                             f"\t{offsets[cid] + pos}\tFracMinHash\t+\t{flags}"
                             f"\t{min(200, round(gc / 50))}\t4294967295\t{cnt}\n")
        print(f"  {sample}: {sum(1 for k in loci if kmc.get(k)):,} observed",
              file=sys.stderr)


if __name__ == "__main__":
    main()
