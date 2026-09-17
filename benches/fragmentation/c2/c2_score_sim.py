#!/usr/bin/env python3
"""C2 Task 2.3 scorer: recall / RMSE / bias / spurious + WCG per arm per cell.

Per cell tag: truth = genomes picked with true_log2ptr (complete-coordinate
simulation). Per arm (frag16 / complete16): a genome is RECALLED when
output.tsv carries a finite log2(PTR) for it; RMSE/bias over recalled
genomes vs true_log2ptr; spurious = estimated genomes not in the truth.
WCG is computed per genome from windows.rates.tsv (fire rate = fraction of
truth genomes flagged). r/slope are absent by design at one sample per cell
(constant-arm rule); the grid-level dose response is read across cov rows.

Writes sim/c2_sim_summary.tsv and prints it.
"""
import glob
import os
import sys

import numpy as np
import pandas as pd

BASE = "/lustre1/g/aos_shihuang/sk2bgrow-hpc"
C2 = os.path.join(BASE, "bench/C2")
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from qc_within_contig import wcg_stat  # noqa: E402

WINDOWS_COLS = ["genome", "enzyme", "contig_id", "start", "end",
                "log2_rate", "log2_se"]


def main():
    rows = []
    for truth_path in sorted(glob.glob(f"{C2}/sim/res/*.truth")):
        tag = os.path.basename(truth_path)[:-6]
        truth = pd.read_csv(truth_path, sep="\t")
        for arm in ("frag16", "complete16"):
            od = f"{C2}/sim/res/{arm}_{tag}"
            opath = os.path.join(od, "output.tsv")
            if not os.path.isfile(opath):
                continue
            out = pd.read_csv(opath, sep="\t")
            est = out.dropna(subset=["log2(PTR)"]).set_index("genome")
            tmap = truth.set_index("genome")["true_log2ptr"]
            both = est.join(tmap, how="inner", rsuffix="_t")
            rec = len(both) / len(truth)
            err = both["log2(PTR)"] - both["true_log2ptr"]
            spur = int((~est.index.isin(truth["genome"])).sum())
            # WCG per truth genome
            n_fire = 0
            wpath = os.path.join(od, "windows.rates.tsv")
            if os.path.isfile(wpath):
                w = pd.read_csv(wpath, sep="\t", usecols=WINDOWS_COLS)
                glen = truth.set_index("genome")["length"]
                for gname, gtruth in tmap.items():
                    gw = w[w["genome"] == gname]
                    if gw.empty or gname not in est.index:
                        continue
                    st = wcg_stat(gw, int(glen[gname]),
                                  float(est.loc[gname, "log2(PTR)"]),
                                  float(est.loc[gname, "se"]))
                    n_fire += int(st["flag"])
            rows.append(dict(
                tag=tag, arm=arm, n_strains=len(truth),
                recall=rec, n_est=len(both), spurious=spur,
                rmse=float(np.sqrt((err ** 2).mean())) if len(err) else np.nan,
                bias=float(err.mean()) if len(err) else np.nan,
                wcg_fire=n_fire, wcg_fire_rate=n_fire / len(truth)))
    df = pd.DataFrame(rows)
    df.to_csv(f"{C2}/sim/c2_sim_summary.tsv", sep="\t", index=False,
              float_format="%.6g")
    pd.set_option("display.width", 200)
    print(df.to_string(index=False))
    print(f"\nwrote {C2}/sim/c2_sim_summary.tsv")


if __name__ == "__main__":
    main()
