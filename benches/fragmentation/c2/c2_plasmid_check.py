#!/usr/bin/env python3
"""C2 Task 1 addendum: why does the jump branch fire on K. pneumoniae?

Fires are confined to the only multi-REPlicon reference (chromosome + 6
plasmids). Hypothesis: chromosome->plasmid boundaries carry a genuine
coverage step (plasmid copy number), which J measures correctly but which
is not coordinate scrambling. Test: recompute J on the 25 K. pneumoniae
cells after dropping every window on a non-primary contig (< 10% of the
largest contig). If the fires vanish, the boundary excess is entirely a
replicon-boundary (copy-number) effect, and the fix direction is "jump
branch on primary-replicon fragments only", not a re-tuned threshold.
"""
import json
import os
import sys

import numpy as np
import pandas as pd

BASE = "/lustre1/g/aos_shihuang/sk2bgrow-hpc"
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from qc_within_contig import boundary_jumps  # noqa: E402

MANIFEST = f"{BASE}/bench/C1b/db/Klebsiella_pneumoniae/manifest.json"
COUNTS = f"{BASE}/bench/C1b/counts"

m = json.load(open(MANIFEST))
clen = {c["id"]: c["length"] for g in m["genomes"] for c in g["contigs"]}
LMAX = max(clen.values())

print(f"contig lengths: {clen}; primary floor = {0.1 * LMAX:,.0f} bp")
rows = []
for cell in sorted(os.listdir(COUNTS)):
    if not cell.startswith(("SRR34095772", "SRR34095773", "SRR34095774",
                            "SRR34095775", "SRR34095776")):
        continue
    wpath = os.path.join(COUNTS, cell, "windows.rates.tsv")
    if not os.path.isfile(wpath):
        continue
    w = pd.read_csv(wpath, sep="\t",
                    usecols=["enzyme", "contig_id", "start", "end",
                             "log2_rate", "log2_se"])
    J_all, _, nb_all = boundary_jumps(w)
    big = w[w["contig_id"].map(clen) >= 0.1 * LMAX]
    J_pri, _, nb_pri = boundary_jumps(big)
    fire_all = bool(np.isfinite(J_all) and J_all > 0.03)
    rows.append((cell, J_all, nb_all, J_pri, nb_pri, fire_all))
    print(f"{cell:22s} J_all={J_all!r:>22} n_b={nb_all:3d}  "
          f"J_primary_only={J_pri!r:>22} n_b={nb_pri:3d}  "
          f"fire_all={fire_all}")
df = pd.DataFrame(rows, columns=["cell", "J_all", "n_boundary_all",
                                 "J_primary", "n_boundary_primary",
                                 "fire_all"])
df.to_csv(f"{BASE}/bench/C2/c2_kpneumo_plasmid_check.tsv", sep="\t",
          index=False, float_format="%.6g")
n_fire = df["fire_all"].sum()
print(f"\nfires with all contigs: {n_fire}/{len(df)}; "
      f"fires primary-contig-only: "
      f"{int((df['J_primary'] > 0.03).sum())}/{len(df)} "
      f"(all J_primary NaN = no within-chromosome boundaries -> quiet)")
