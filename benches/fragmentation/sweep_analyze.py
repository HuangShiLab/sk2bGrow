#!/usr/bin/env python3
"""Where does the coordinate stop being usable? Accuracy against contig count.

    WORK=/path/to/work python3 sweep_analyze.py
"""
import glob
import os
import re
import sys

import numpy as np
import pandas as pd
from scipy import stats


def _n50(layout):
    """N50 of the contig set, read back from fragment.py's truth file."""
    if not os.path.exists(layout):
        return float("nan")
    lens = np.sort(pd.read_csv(layout, sep="\t")["length"].to_numpy())[::-1]
    return float(lens[np.searchsorted(np.cumsum(lens), lens.sum() / 2)])


HERE = os.path.dirname(os.path.abspath(__file__))
WORK = os.environ.get("WORK") or sys.exit("set WORK")
gt = pd.read_csv(os.path.join(HERE, "..", "zheng2020", "growth_rates.tsv"),
                 sep="\t").set_index("medium")

rows = []
for d in sorted(glob.glob(os.path.join(WORK, "sweep", "out", "*/"))):
    m = re.match(r"^n(\d+)_(.+?)_([\d.]+)x$", os.path.basename(d.rstrip("/")))
    if not m:
        continue
    f = os.path.join(d, "output.tsv")
    if not os.path.exists(f):
        continue
    df = pd.read_csv(f, sep="\t", na_values=["NA", "n/a"])
    if df.empty:
        continue
    rows.append(dict(n_contigs=int(m.group(1)), medium=m.group(2),
                     cov=float(m.group(3)), log2ptr=df.iloc[0]["log2(PTR)"],
                     passed=bool(df.iloc[0].get("pass_qc", True))))

res = pd.DataFrame(rows)
if res.empty:
    sys.exit("no sweep results")
res["growth_rate"] = res["medium"].map(gt["growth_rate"])
res["pred_log2ptr"] = res["medium"].map(gt["pred_log2ptr"])
res.to_csv(os.path.join(WORK, "sweep_raw.tsv"), sep="\t", index=False,
           float_format="%.4f")

grow = res[res["medium"] != "RUN_OUT"]
cov = grow["cov"].iloc[0]
print(f"CONTIG-COUNT SWEEP at {cov:g}x  (E. coli, 4.64 Mbp, {grow['medium'].nunique()} media)")
print(f"{'contigs':>8}{'N50 (kb)':>10}{'n':>4}{'Pearson r':>11}{'RMSE':>8}"
      f"{'bias':>8}{'slope':>8}{'QC pass':>9}")
for n in sorted(grow["n_contigs"].unique()):
    s = grow[grow["n_contigs"] == n].dropna(subset=["log2ptr"])
    if len(s) < 3:
        continue
    ok = s.dropna(subset=["pred_log2ptr"])
    e = ok["log2ptr"] - ok["pred_log2ptr"]
    n50 = _n50(os.path.join(WORK, "sweep", f"n{n}.layout.tsv"))
    print(f"{n:>8}{n50 / 1000:>10,.0f}{len(s):>4}"
          f"{stats.pearsonr(s['growth_rate'], s['log2ptr'])[0]:>11.4f}"
          f"{np.sqrt((e ** 2).mean()):>8.3f}{e.mean():>8.3f}"
          f"{np.polyfit(s['growth_rate'], s['log2ptr'], 1)[0]:>8.3f}"
          f"{100 * s['passed'].mean():>8.0f}%")
