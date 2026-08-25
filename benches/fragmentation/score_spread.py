#!/usr/bin/env python3
"""Score the order-free spread estimator against the V-fit and Pilea.

Reads the windows.rates.tsv each profile run already wrote, so it needs no
re-run of the pipeline.

    WORK=/path/to/work python3 score_spread.py
"""
import glob
import os
import re
import sys

import numpy as np
import pandas as pd
from scipy import stats

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from spread_estimator import fuse, per_enzyme

HERE = os.path.dirname(os.path.abspath(__file__))
WORK = os.environ.get("WORK") or sys.exit("set WORK")
gt = pd.read_csv(os.path.join(HERE, "..", "zheng2020", "growth_rates.tsv"),
                 sep="\t").set_index("medium")

rows = []
for d in sorted(glob.glob(os.path.join(WORK, "out", "*/"))):
    m = re.match(r"^(complete|frag|scafSelf|scafRel)_(.+?)_([\d.]+)x$",
                 os.path.basename(d.rstrip("/")))
    f = os.path.join(d, "windows.rates.tsv")
    if not m or not os.path.exists(f):
        continue
    est, se, Q, k = fuse(per_enzyme(pd.read_csv(f, sep="\t")))
    rows.append(dict(cond="spread_" + m.group(1), medium=m.group(2),
                     cov=float(m.group(3)), log2ptr=est, spread_se=se,
                     Q=Q, n_enzymes=k))
sp = pd.DataFrame(rows)
sp.to_csv(os.path.join(WORK, "spread_raw.tsv"), sep="\t", index=False,
          float_format="%.4f")

vf = pd.read_csv(os.path.join(WORK, "results_raw.tsv"), sep="\t")
d = pd.concat([vf[["cond", "medium", "cov", "log2ptr"]],
               sp[["cond", "medium", "cov", "log2ptr"]]])
d["growth_rate"] = d["medium"].map(gt["growth_rate"])
d["pred"] = d["medium"].map(gt["pred_log2ptr"])

LAB = {"complete": "V-fit, complete (reference)",
       "frag": "V-fit on 100 contigs",
       "spread_frag": "spread-MLE on 100 contigs",
       "pileaFrag": "Pilea on 100 contigs"}
g = d[d["medium"] != "RUN_OUT"].dropna(subset=["log2ptr"])
print("ON A FRAGMENTED REFERENCE — which estimator survives?")
print(f"{'cov':>5}  {'estimator':28}{'n':>4}{'r':>8}{'RMSE':>8}{'bias':>8}{'slope':>8}")
for cov in sorted(g["cov"].unique()):
    for c in LAB:
        s = g[(g["cov"] == cov) & (g["cond"] == c)]
        if len(s) < 3 or s["log2ptr"].nunique() < 2:
            print(f"{cov:>5g}  {LAB[c]:28}{len(s):>4}{'--':>8}{'--':>8}{'--':>8}{'--':>8}")
            continue
        e = s["log2ptr"] - s["pred"]
        print(f"{cov:>5g}  {LAB[c]:28}{len(s):>4}"
              f"{stats.pearsonr(s['growth_rate'], s['log2ptr'])[0]:>8.3f}"
              f"{np.sqrt((e ** 2).mean()):>8.3f}{e.mean():>8.3f}"
              f"{np.polyfit(s['growth_rate'], s['log2ptr'], 1)[0]:>8.3f}")
    print()

print("NEGATIVE CONTROL — RUN_OUT, stationary phase, truth ~ 0")
ctl = d[(d["medium"] == "RUN_OUT") & d["cond"].isin(LAB)]
print(ctl.pivot_table(index="cov", columns="cond", values="log2ptr")
      .round(3).to_string())
